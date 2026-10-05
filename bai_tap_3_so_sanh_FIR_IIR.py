# -*- coding: utf-8 -*-
"""Bài tập 3: So sánh MA (FIR) và EMA (IIR) bằng mô phỏng Monte Carlo."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_SNR_DB = np.arange(-10, 21, 1)


@dataclass(frozen=True)
class MoPhongConfig:
    N: int = 600
    f0: float = 0.01
    M: int = 5
    alpha: float = 0.2
    snr_start: int = -10
    snr_stop: int = 20
    snr_step: int = 1
    num_trials: int = 100
    seed: int = 42
    n_transient: int = 50
    ema_init: str = "zero"  # zero | first_sample
    output_dir: Path = Path("outputs/bai_tap_3")

    def snr_range(self) -> np.ndarray:
        return np.arange(self.snr_start, self.snr_stop + 1, self.snr_step)


def tao_tin_hieu_sach(N: int, f0: float) -> Tuple[np.ndarray, np.ndarray]:
    n = np.arange(N)
    s = np.sin(2 * np.pi * f0 * n)
    return n, s


def sinh_nhieu_awgn(signal: np.ndarray, snr_db: float, rng: np.random.Generator) -> np.ndarray:
    p_signal = float(np.mean(signal**2))
    p_noise = p_signal / (10.0 ** (snr_db / 10.0))
    sigma = np.sqrt(p_noise)
    return rng.normal(0.0, sigma, size=signal.shape)


def loc_ma(x: np.ndarray, M: int) -> np.ndarray:
    kernel = np.ones(M, dtype=float) / float(M)
    y_full = np.convolve(x, kernel, mode="full")
    return y_full[: x.size]


def loc_ema(x: np.ndarray, alpha: float, init_mode: str = "zero") -> np.ndarray:
    y = np.empty_like(x, dtype=float)
    y_prev = 0.0 if init_mode == "zero" else float(x[0])
    for i, xi in enumerate(x):
        y_prev = alpha * xi + (1.0 - alpha) * y_prev
        y[i] = y_prev
    return y


def tinh_rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def do_tre_ma(M: int) -> float:
    return (M - 1) / 2.0


def do_tre_ema_xap_xi_dc(alpha: float) -> float:
    return (1.0 - alpha) / alpha


def cat_vung_cong_bang(
    y: np.ndarray,
    s: np.ndarray,
    delay_samples: float,
    n_transient: int,
    max_delay_samples: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cắt cùng vùng chỉ số sạch cho mọi bộ lọc.

    So sánh trên chỉ số sạch n thuộc [n_transient, N-max_delay-1].
    Với bộ lọc có trễ d, tín hiệu so sánh là y[n + d] so với s[n].
    """
    d = int(round(delay_samples))
    max_d = int(max_delay_samples)
    n0 = n_transient
    n1 = s.size - max_d
    if n1 <= n0:
        raise ValueError("Cấu hình n_transient/max_delay làm rỗng vùng đánh giá.")
    idx_s = np.arange(n0, n1)
    idx_y = idx_s + d
    return y[idx_y], s[idx_s], idx_s


def thong_ke_monte_carlo(cfg: MoPhongConfig) -> Dict[str, np.ndarray]:
    rng = np.random.default_rng(cfg.seed)
    _, s_clean = tao_tin_hieu_sach(cfg.N, cfg.f0)

    snr_db = cfg.snr_range()
    delay_ma = do_tre_ma(cfg.M)
    delay_ema = do_tre_ema_xap_xi_dc(cfg.alpha)
    max_delay = int(np.ceil(max(delay_ma, delay_ema)))

    num_snr = snr_db.size
    rmse_raw = np.zeros((num_snr, cfg.num_trials), dtype=float)
    rmse_ma_aligned = np.zeros_like(rmse_raw)
    rmse_ema_aligned = np.zeros_like(rmse_raw)
    rmse_ma_unaligned = np.zeros_like(rmse_raw)
    rmse_ema_unaligned = np.zeros_like(rmse_raw)

    idx_common = np.arange(cfg.n_transient, cfg.N - max_delay)

    for i, snr in enumerate(snr_db):
        for t in range(cfg.num_trials):
            noise = sinh_nhieu_awgn(s_clean, float(snr), rng)
            x_noisy = s_clean + noise
            y_ma = loc_ma(x_noisy, cfg.M)
            y_ema = loc_ema(x_noisy, cfg.alpha, init_mode=cfg.ema_init)

            # Metric chính: RMSE sau căn chỉnh trễ + bỏ quá độ + cùng vùng chỉ số
            y_ma_a, s_ma_a, _ = cat_vung_cong_bang(y_ma, s_clean, delay_ma, cfg.n_transient, max_delay)
            y_ema_a, s_ema_a, _ = cat_vung_cong_bang(y_ema, s_clean, delay_ema, cfg.n_transient, max_delay)
            y_raw_a, s_raw_a, _ = cat_vung_cong_bang(x_noisy, s_clean, 0.0, cfg.n_transient, max_delay)

            rmse_ma_aligned[i, t] = tinh_rmse(y_ma_a, s_ma_a)
            rmse_ema_aligned[i, t] = tinh_rmse(y_ema_a, s_ema_a)
            rmse_raw[i, t] = tinh_rmse(y_raw_a, s_raw_a)

            # Metric phụ: chưa căn chỉnh để minh họa ảnh hưởng độ trễ
            rmse_ma_unaligned[i, t] = tinh_rmse(y_ma[idx_common], s_clean[idx_common])
            rmse_ema_unaligned[i, t] = tinh_rmse(y_ema[idx_common], s_clean[idx_common])

    return {
        "snr_db": snr_db,
        "rmse_raw_mean": rmse_raw.mean(axis=1),
        "rmse_raw_std": rmse_raw.std(axis=1),
        "rmse_ma_aligned_mean": rmse_ma_aligned.mean(axis=1),
        "rmse_ma_aligned_std": rmse_ma_aligned.std(axis=1),
        "rmse_ema_aligned_mean": rmse_ema_aligned.mean(axis=1),
        "rmse_ema_aligned_std": rmse_ema_aligned.std(axis=1),
        "rmse_ma_unaligned_mean": rmse_ma_unaligned.mean(axis=1),
        "rmse_ema_unaligned_mean": rmse_ema_unaligned.mean(axis=1),
        "delay_ma": np.array([delay_ma], dtype=float),
        "delay_ema_approx_lowfreq": np.array([delay_ema], dtype=float),
        "max_delay_eval": np.array([max_delay], dtype=float),
    }


def phan_tich_tu_du_lieu(res: Dict[str, np.ndarray]) -> Dict[str, object]:
    diff = res["rmse_ema_aligned_mean"] - res["rmse_ma_aligned_mean"]
    snr = res["snr_db"]
    ema_tot_hon = snr[diff < 0]
    ma_tot_hon = snr[diff > 0]
    hoa_nhau = snr[np.isclose(diff, 0.0)]

    cross_idx = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
    crossing = []
    for i in cross_idx:
        crossing.append((int(snr[i]), int(snr[i + 1])))

    return {
        "ema_better_snr": ema_tot_hon.tolist(),
        "ma_better_snr": ma_tot_hon.tolist(),
        "tie_snr": hoa_nhau.tolist(),
        "crossings": crossing,
    }


def ve_hinh_1_input_filter_output(cfg: MoPhongConfig, out_dir: Path) -> None:
    rng = np.random.default_rng(cfg.seed)
    n, s_clean = tao_tin_hieu_sach(N=220, f0=0.015)
    snr_demo = 3
    x_noisy = s_clean + sinh_nhieu_awgn(s_clean, snr_demo, rng)
    y_ma = loc_ma(x_noisy, cfg.M)
    y_ema = loc_ema(x_noisy, cfg.alpha, cfg.ema_init)

    fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
    fig.suptitle("(a) Input → Filter → Output", fontweight="bold")

    axes[0].plot(n, s_clean, "k--", label="s[n] sạch")
    axes[0].plot(n, x_noisy, color="#e67e22", alpha=0.7, label=f"x[n] (AWGN, SNR={snr_demo} dB)")
    axes[0].legend(loc="upper right")
    axes[0].grid(ls=":", alpha=0.5)

    axes[1].plot(n, s_clean, "k--", label="s[n]")
    axes[1].plot(n, y_ma, "b", label=f"MA (M={cfg.M})")
    axes[1].plot(n, y_ema, "r", label=f"EMA (α={cfg.alpha})")
    axes[1].legend(loc="upper right")
    axes[1].grid(ls=":", alpha=0.5)

    axes[2].plot(n[30:110], s_clean[30:110], "k", lw=2, label="s[n]")
    axes[2].plot(n[30:110], y_ma[30:110], "b", label=f"MA (trễ {do_tre_ma(cfg.M):.1f})")
    axes[2].plot(n[30:110], y_ema[30:110], "r", label=f"EMA (trễ thấp tần xấp xỉ {do_tre_ema_xap_xi_dc(cfg.alpha):.1f})")
    axes[2].set_xlabel("n")
    axes[2].legend(loc="lower left")
    axes[2].grid(ls=":", alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_dir / "slide_hinh1_input_filter_output.png", dpi=160)
    plt.close(fig)


def ve_hinh_2_nguyen_ly(cfg: MoPhongConfig, out_dir: Path) -> None:
    n = np.arange(25)
    h_ma = np.zeros_like(n, dtype=float)
    h_ma[: cfg.M] = 1.0 / cfg.M
    h_ema = cfg.alpha * (1.0 - cfg.alpha) ** n

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].stem(n, h_ma, basefmt="k-")
    axes[0].set_title(f"(b) MA/FIR, M={cfg.M}")
    axes[0].grid(ls=":", alpha=0.5)

    axes[1].stem(n, h_ema, basefmt="k-")
    axes[1].set_title(f"(c) EMA/IIR, α={cfg.alpha}")
    axes[1].grid(ls=":", alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_dir / "slide_hinh2_nguyen_ly_dap_ung_xung.png", dpi=160)
    plt.close(fig)


def ve_hinh_3_rmse(res: Dict[str, np.ndarray], cfg: MoPhongConfig, out_dir: Path) -> None:
    snr = res["snr_db"]
    ma = res["rmse_ma_aligned_mean"]
    ema = res["rmse_ema_aligned_mean"]
    raw = res["rmse_raw_mean"]
    diff_pct = (ema - ma) / ma * 100.0

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8), gridspec_kw={"height_ratios": [2.2, 1]})
    ax1.plot(snr, raw, "k:", label="Raw")
    ax1.plot(snr, ma, "b-o", ms=4, label="MA aligned")
    ax1.plot(snr, ema, "r-s", ms=4, label="EMA aligned")
    ax1.set_ylabel("RMSE")
    ax1.set_title("(e) RMSE theo SNR (metric chính: aligned RMSE)")
    ax1.grid(ls=":", alpha=0.5)
    ax1.legend()

    min_snr_idx = int(np.argmin(snr))
    max_snr_idx = int(np.argmax(snr))
    winner_low = "EMA" if ema[min_snr_idx] < ma[min_snr_idx] else "MA"
    winner_high = "EMA" if ema[max_snr_idx] < ma[max_snr_idx] else "MA"
    ax1.text(
        0.01,
        0.02,
        f"SNR thấp nhất ({int(snr[min_snr_idx])} dB): {winner_low} tốt hơn\n"
        f"SNR cao nhất ({int(snr[max_snr_idx])} dB): {winner_high} tốt hơn",
        transform=ax1.transAxes,
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    ax2.bar(snr, diff_pct, color=["#c0392b" if v < 0 else "#2980b9" for v in diff_pct])
    ax2.axhline(0, color="k", ls="--", lw=1)
    ax2.set_xlabel("SNR (dB)")
    ax2.set_ylabel("(EMA-MA)/MA (%)")
    ax2.grid(ls=":", alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_dir / "slide_hinh3_so_sanh_rmse_snr.png", dpi=160)
    plt.close(fig)


def ve_hinh_4_tradeoff(cfg: MoPhongConfig, out_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    var_ma = 1.0 / cfg.M
    var_ema = cfg.alpha / (2.0 - cfg.alpha)
    axes[0, 0].bar(["MA", "EMA"], [var_ma, var_ema], color=["#2980b9", "#c0392b"])
    axes[0, 0].set_title("Giảm phương sai nhiễu trắng (noise-only)")
    axes[0, 0].set_ylabel("σ²_out/σ²_in")
    axes[0, 0].grid(ls=":", alpha=0.5)

    delay_ma = do_tre_ma(cfg.M)
    delay_ema = do_tre_ema_xap_xi_dc(cfg.alpha)
    axes[0, 1].bar(["MA", "EMA"], [delay_ma, delay_ema], color=["#2980b9", "#c0392b"])
    axes[0, 1].set_title("Độ trễ: MA cố định, EMA xấp xỉ thấp tần")
    axes[0, 1].set_ylabel("Mẫu")
    axes[0, 1].grid(ls=":", alpha=0.5)

    w0 = 2.0 * np.pi * cfg.f0
    mag_ma = np.abs(np.sin(cfg.M * w0 / 2.0) / (cfg.M * np.sin(w0 / 2.0)))
    mag_ema = cfg.alpha / np.abs(1.0 - (1.0 - cfg.alpha) * np.exp(-1j * w0))
    att_ma = (1.0 - mag_ma) * 100.0
    att_ema = (1.0 - mag_ema) * 100.0
    axes[1, 0].bar(["MA", "EMA"], [att_ma, att_ema], color=["#2980b9", "#c0392b"])
    axes[1, 0].set_title(f"Over-smoothing tại f0={cfg.f0}")
    axes[1, 0].set_ylabel("Suy hao biên độ (%)")
    axes[1, 0].grid(ls=":", alpha=0.5)

    f = np.linspace(1e-4, 0.15, 300)
    w = 2 * np.pi * f
    beta = 1.0 - cfg.alpha
    tau_ma = np.full_like(f, delay_ma)
    tau_ema = beta * (np.cos(w) - beta) / (1 + beta**2 - 2 * beta * np.cos(w))
    axes[1, 1].plot(f, tau_ma, "b", label="MA")
    axes[1, 1].plot(f, tau_ema, "r", label="EMA")
    axes[1, 1].axvline(cfg.f0, color="green", ls=":", label="f0")
    axes[1, 1].set_title("Độ trễ nhóm theo tần số")
    axes[1, 1].set_xlabel("f (cycles/sample)")
    axes[1, 1].set_ylabel("Mẫu")
    axes[1, 1].grid(ls=":", alpha=0.5)
    axes[1, 1].legend()

    fig.tight_layout()
    fig.savefig(out_dir / "slide_hinh4_phan_tich_danh_doi_va_han_che.png", dpi=160)
    plt.close(fig)


def xuat_csv_json(cfg: MoPhongConfig, res: Dict[str, np.ndarray], analysis: Dict[str, object], out_dir: Path) -> None:
    rows = []
    for i, snr in enumerate(res["snr_db"]):
        ma = float(res["rmse_ma_aligned_mean"][i])
        ema = float(res["rmse_ema_aligned_mean"][i])
        rows.append(
            {
                "snr_db": int(snr),
                "rmse_raw_mean": float(res["rmse_raw_mean"][i]),
                "rmse_ma_aligned_mean": ma,
                "rmse_ema_aligned_mean": ema,
                "rmse_ma_unaligned_mean": float(res["rmse_ma_unaligned_mean"][i]),
                "rmse_ema_unaligned_mean": float(res["rmse_ema_unaligned_mean"][i]),
                "rmse_ma_aligned_std": float(res["rmse_ma_aligned_std"][i]),
                "rmse_ema_aligned_std": float(res["rmse_ema_aligned_std"][i]),
                "bo_loc_tot_hon": "EMA" if ema < ma else "MA",
                "chenh_lech_pct_(EMA-MA)/MA": float((ema - ma) / ma * 100.0),
            }
        )

    csv_path = out_dir / "ket_qua_rmse_theo_snr.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    json_path = out_dir / "tong_hop_ket_qua.json"
    payload = {
        "config": {**asdict(cfg), "output_dir": str(cfg.output_dir)},
        "analysis": analysis,
        "results": rows,
        "ly_thuyet": {
            "ma_noise_variance_ratio": 1.0 / cfg.M,
            "ema_noise_variance_ratio": cfg.alpha / (2.0 - cfg.alpha),
            "ghi_chu": "Các hệ số trên chỉ áp dụng cho thành phần nhiễu trắng; RMSE tổng còn chịu ảnh hưởng bias và đáp ứng biên độ.",
        },
    }
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def in_tom_tat_console(cfg: MoPhongConfig, res: Dict[str, np.ndarray], analysis: Dict[str, object]) -> None:
    print("\n=== TÓM TẮT KẾT QUẢ ===")
    print(f"M={cfg.M}, alpha={cfg.alpha}, trials={cfg.num_trials}, seed={cfg.seed}")
    print(
        f"Delay MA={res['delay_ma'][0]:.1f} mẫu (chính xác tuyến tính), "
        f"Delay EMA thấp tần≈{res['delay_ema_approx_lowfreq'][0]:.1f} mẫu (xấp xỉ, phụ thuộc tần số)."
    )
    if analysis["crossings"]:
        print(f"Các khoảng giao cắt MA/EMA theo SNR: {analysis['crossings']}")
    else:
        print("Không có giao cắt MA/EMA trong dải SNR đã quét.")


def chay_self_check() -> None:
    cfg = MoPhongConfig(num_trials=5, N=300)
    rng = np.random.default_rng(cfg.seed)
    _, s = tao_tin_hieu_sach(cfg.N, cfg.f0)

    # 1) MA coefficients sum = 1
    kernel = np.ones(cfg.M) / cfg.M
    assert np.isclose(kernel.sum(), 1.0), "Tổng hệ số MA phải bằng 1"

    # 2) EMA impulse response gần 1
    n = np.arange(3000)
    h_ema = cfg.alpha * (1 - cfg.alpha) ** n
    assert np.isclose(h_ema.sum(), 1.0, atol=1e-6), "Tổng đáp ứng xung EMA phải xấp xỉ 1"

    # 3) SNR tạo ra gần giá trị yêu cầu
    target_snr = 5
    noise = sinh_nhieu_awgn(s, target_snr, rng)
    snr_meas = 10 * np.log10(np.mean(s**2) / np.mean(noise**2))
    assert abs(snr_meas - target_snr) < 1.0, "SNR sinh ra lệch quá mức"

    # 4) shape đầu ra
    x = s + noise
    y_ma = loc_ma(x, cfg.M)
    y_ema = loc_ema(x, cfg.alpha)
    assert y_ma.shape == x.shape == y_ema.shape, "Shape output phải khớp input"

    print("Self-check: PASS")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="So sánh MA/EMA theo RMSE-SNR (Monte Carlo).")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/bai_tap_3"))
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n", type=int, default=600)
    parser.add_argument("--f0", type=float, default=0.01)
    parser.add_argument("--m", type=int, default=5)
    parser.add_argument("--alpha", type=float, default=0.2)
    parser.add_argument("--snr-start", type=int, default=-10)
    parser.add_argument("--snr-stop", type=int, default=20)
    parser.add_argument("--snr-step", type=int, default=1)
    parser.add_argument("--transient", type=int, default=50)
    parser.add_argument("--ema-init", choices=["zero", "first_sample"], default="zero")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--self-check", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_check:
        chay_self_check()
        return

    cfg = MoPhongConfig(
        N=args.n,
        f0=args.f0,
        M=args.m,
        alpha=args.alpha,
        snr_start=args.snr_start,
        snr_stop=args.snr_stop,
        snr_step=args.snr_step,
        num_trials=args.trials,
        seed=args.seed,
        n_transient=args.transient,
        ema_init=args.ema_init,
        output_dir=args.output_dir,
    )

    out_dir = cfg.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    res = thong_ke_monte_carlo(cfg)
    analysis = phan_tich_tu_du_lieu(res)
    xuat_csv_json(cfg, res, analysis, out_dir)

    if not args.skip_figures:
        ve_hinh_1_input_filter_output(cfg, out_dir)
        ve_hinh_2_nguyen_ly(cfg, out_dir)
        ve_hinh_3_rmse(res, cfg, out_dir)
        ve_hinh_4_tradeoff(cfg, out_dir)

    in_tom_tat_console(cfg, res, analysis)
    print(f"Đã lưu kết quả tại: {out_dir.resolve()}")


if __name__ == "__main__":
    main()

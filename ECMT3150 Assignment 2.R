
# ECMT3150 Assignment 2 (Semester 2, 2025) -- R script
# Numerical answers are rounded to 3 decimal places, as required.


r3 <- function(x) round(x, 3)
cat_line <- function(label, x) cat(sprintf("%-55s %s\n", label, format(x, nsmall = 3)))


# Question 1: Single-period binomial model for NVDA

S0 <- 100; Su <- 120; Sd <- 80; r <- 0.1; K <- 100

cat("\n=== Question 1 ===\n")


# (a) Martingale condition: S0 = (1/(1+r)) * [p*Su + (1-p)*Sd]
#     => p = ((1+r)*S0 - Sd) / (Su - Sd)
p <- ((1 + r) * S0 - Sd) / (Su - Sd)
cat_line("(a) Risk-neutral probability p:", r3(p))
stopifnot(abs(p - 0.75) < 1e-12)


# (b) European call price: discounted risk-neutral expectation of payoff
Cu <- max(Su - K, 0); Cd <- max(Sd - K, 0)
C0 <- (p * Cu + (1 - p) * Cd) / (1 + r)
cat_line("(b) Fair call price C0:", r3(C0))


# (c) Replicating portfolio: D units of money account (worth 1 at t=0, 1+r at t=1),
#     Delta shares of NVDA.
#       D*(1+r) + Delta*Su = Cu
#       D*(1+r) + Delta*Sd = Cd
Delta <- (Cu - Cd) / (Su - Sd)
D     <- (Cu - Delta * Su) / (1 + r)
cat_line("(c) Delta (shares of NVDA, + = long):", r3(Delta))
cat_line("    D (units of money account, - = short/borrow):", r3(D))
cat_line("    Portfolio cost D + Delta*S0 (should equal C0):", r3(D + Delta * S0))
stopifnot(abs((D + Delta * S0) - C0) < 1e-10)


# (d)(i) Market price of the call is $15 > fair price: call is overpriced.
#        Strategy: sell n calls, buy replicating portfolio per call, pocket difference.
Cmkt   <- 15
profit_per_call <- Cmkt - C0
n_calls <- 300000 / profit_per_call
cat_line("(d)(i) Mispricing per call (Cmkt - C0):", r3(profit_per_call))
cat_line("      Number of calls to sell for $300,000:", r3(n_calls))
cat_line("      Shares of NVDA to buy (n * Delta):", r3(n_calls * Delta))
cat_line("      Money-account position (n * D, - = borrow):", r3(n_calls * D))

# Payoff check at time 1 in both states: the position nets to exactly 0 time-1 cash
# and the $300,000 is earned upfront at t=0 (financed at risk-free rate).
payoff_check <- function(ST) {
  n_calls * (-max(ST - K, 0) + Delta * ST + D * (1 + r))
}
cat_line("      Time-1 net payoff if up   (should be 0):", r3(payoff_check(Su)))
cat_line("      Time-1 net payoff if down (should be 0):", r3(payoff_check(Sd)))
cat_line("      Time-0 cash inflow:", r3(n_calls * (Cmkt - (D + Delta * S0))))


# (d)(ii) Discussion points (see written answer): e.g. transaction costs / bid-ask
#         spreads, borrowing/short-selling constraints and margin requirements,
#         position/liquidity limits, discrete trading and model risk.


# Question 2: AUD/CAD one-period model (r = 0) and diffusion model

cat("\n=== Question 2 ===\n")


# (a) E_Q[S1] = S0  =>  q*S0*u + (1-q)*S0/u = S0  =>  q = 1/(u+1)
q_fun <- function(u) 1 / (u + 1)
# Numeric verification of the martingale property for a few u values
for (u in c(1.05, 1.2, 1.5)) {
  S0b <- 1
  Eq  <- q_fun(u) * S0b * u + (1 - q_fun(u)) * S0b / u
  cat(sprintf("(a) u = %.2f: q = %.3f, E_Q[S1] = %.3f\n", u, q_fun(u), Eq))
}


# (b) ATM call on AUD/CAD, S0 = 1, K = 1, r = 0:  payoff (u - 1) w.p. q, else 0
call_AUDCAD <- function(u) q_fun(u) * (u - 1)          # = (u-1)/(u+1)


# (c) ATM call on CAD/AUD, S0 = 1, K = 1: CAD/AUD = 1/S_t.
#     Up-move of S (prob q)   : 1/S1 = 1/u  -> payoff 0
#     Down-move of S (1-q)    : 1/S1 = u    -> payoff (u - 1)
#     Following the hint (same Q), price = (1-q)*(u-1) = u(u-1)/(u+1).
call_CADAUD <- function(u) (1 - q_fun(u)) * (u - 1)

# Remark: u(u-1)/(u+1) = u * [(u-1)/(u+1)], i.e. the CAD/AUD call is worth u times
# the AUD/CAD call under the same measure Q (since 1/S_t is NOT a Q-martingale;
# this is the Siegel/Jensen effect). Under the CAD/AUD-numeraire measure Q' the
# two prices would coincide at (u-1)/(u+1) by symmetry.
u_grid <- c(1.05, 1.1, 1.2, 1.5, 2)
tab <- data.frame(
  u            = u_grid,
  q            = r3(q_fun(u_grid)),
  call_AUDCAD  = r3(call_AUDCAD(u_grid)),
  call_CADAUD  = r3(call_CADAUD(u_grid)),
  ratio        = r3(call_CADAUD(u_grid) / call_AUDCAD(u_grid))
)
cat("\n(b)-(c) Prices at selected u (symbolic: (u-1)/(u+1) and u(u-1)/(u+1)):\n")
print(tab, row.names = FALSE)


# (d) Ito's lemma. dS = sigma S dW (zero drift under Q), Y = f(S) = S^{-1}
#     f' = -S^{-2}, f'' = 2 S^{-3}, no explicit time dependence:
#     dY = f' dS + 0.5 f'' (dS)^2
#        = -S^{-2} (sigma S dW) + 0.5 * 2 S^{-3} * sigma^2 S^2 dt
#        = sigma^2 Y dt - sigma Y dW
# Symbolic check with base R's deriv()
f <- expression(S^(-1))
cat("\n(d) Symbolic derivatives of Y = 1/S via deriv():\n")
d1 <- D(f[[1]], "S");  d2 <- D(d1, "S")
cat("    dY/dS     =", deparse(d1), "\n")
cat("    d2Y/dS2   =", deparse(d2), "\n")
cat("    => dY = sigma^2 * Y dt - sigma * Y dW\n")

# Monte Carlo check: under Q, S(t) = S0*exp(-sigma^2 t/2 + sigma W_t).
# From the derived SDE, E[Y(t)] = Y(0) * exp(sigma^2 t).
set.seed(3150)
sigma <- 0.2; Tt <- 1; nsim <- 1e6
Wt <- rnorm(nsim, 0, sqrt(Tt))
St <- 1 * exp(-0.5 * sigma^2 * Tt + sigma * Wt)
Yt <- 1 / St
cat_line("    MC  E[Y(T)] (sigma = 0.2, T = 1):", r3(mean(Yt)))
cat_line("    Theory exp(sigma^2 T):", r3(exp(sigma^2 * Tt)))
cat_line("    MC  E[S(T)] (should be 1, martingale):", r3(mean(St)))


# (e) Flash crash: S(t) -> 0  =>  Y = 1/S -> infinity.
#     Drift = sigma^2 * Y(t) -> +inf; absolute volatility = sigma * Y(t) -> +inf
#     (the volatility is constant only in proportional terms: sigma).
cat("\n(e) Drift and absolute volatility of Y as S falls:\n")
S_path <- c(1, 0.5, 0.1, 0.01, 0.001)
Ydf <- data.frame(
  S         = S_path,
  Y         = 1 / S_path,
  drift     = r3(sigma^2 / S_path),
  abs_vol   = r3(sigma / S_path)
)
print(Ydf, row.names = FALSE)


# FIGURES (base R only -- no extra packages needed)
# Each figure is saved as a 300-dpi PNG and a vector PDF in ./figures/

cat("\n=== Figures ===\n")

fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

pal <- list(blue = "#1F4E79", red = "#C0392B", green = "#2E7D32",
            orange = "#E67E22", grey = "#7F8C8D", lightgrey = "#D5D8DC")

# Draw one figure into both a PNG and a PDF device.
save_fig <- function(name, plot_fun, width = 9, height = 5.5) {
  png(file.path(fig_dir, paste0(name, ".png")),
      width = width, height = height, units = "in", res = 300)
  plot_fun(); invisible(dev.off())
  pdf(file.path(fig_dir, paste0(name, ".pdf")), width = width, height = height)
  plot_fun(); invisible(dev.off())
  cat("  saved:", file.path(fig_dir, name), "(.png/.pdf)\n")
}


# ---- Figure 1: Q1 binomial tree ----------------------------------------------
fig_tree <- function() {
  par(mar = c(1, 1, 3, 1))
  plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10))
  title(main = "Q1: One-period binomial model for NVDA (risk-neutral measure Q)")
  node <- function(x, y, txt, fill) {
    rect(x - 1.7, y - 1.2, x + 1.7, y + 1.2, col = fill, border = pal$blue, lwd = 2)
    text(x, y, txt, cex = 1)
  }
  node(1.9, 5, sprintf("t = 0\nS0 = %.0f\nC0 = %.3f", S0, C0), "#EAF2F8")
  node(8.1, 8, sprintf("t = 1 (up)\nS1 = %.0f\nCall payoff = %.0f", Su, Cu), "#E8F5E9")
  node(8.1, 2, sprintf("t = 1 (down)\nS1 = %.0f\nCall payoff = %.0f", Sd, Cd), "#FDEDEC")
  arrows(3.7, 5.6, 6.3, 7.4, length = 0.15, lwd = 2, col = pal$green)
  arrows(3.7, 4.4, 6.3, 2.6, length = 0.15, lwd = 2, col = pal$red)
  text(4.4, 7.0, sprintf("p = %.3f", p), col = pal$green, font = 2)
  text(4.4, 3.0, sprintf("1 - p = %.3f", 1 - p), col = pal$red, font = 2)
  text(5, 0.5, sprintf(
    "Replicating portfolio: Delta = %.3f shares, D = %.3f in the money account (r = %.1f)",
    Delta, D, r), cex = 0.9)
}
save_fig("fig1_binomial_tree", fig_tree, 9, 5.5)


# ---- Figure 2: Q1 replication payoff ------------------------------------------
fig_replication <- function() {
  par(mar = c(4.5, 4.5, 3, 1))
  ST <- seq(60, 140, length.out = 200)
  call_pay <- pmax(ST - K, 0)
  repl_pay <- Delta * ST + D * (1 + r)
  plot(ST, call_pay, type = "l", lwd = 3, col = pal$blue,
       xlab = expression(S[1] ~ "(NVDA price at time 1)"),
       ylab = "Payoff at time 1 ($)",
       main = "Q1(c): Call payoff vs replicating portfolio",
       ylim = range(c(call_pay, repl_pay)))
  lines(ST, repl_pay, lwd = 3, lty = 2, col = pal$orange)
  abline(h = 0, col = pal$lightgrey)
  points(c(Sd, Su), c(Cd, Cu), pch = 19, cex = 1.8, col = pal$red)
  text(Sd, Cd, "down node", pos = 3, col = pal$red, cex = 0.9)
  text(Su, Cu, "up node", pos = 2, col = pal$red, cex = 0.9)
  legend("topleft", bty = "n", lwd = c(3, 3, NA), lty = c(1, 2, NA), pch = c(NA, NA, 19),
         col = c(pal$blue, pal$orange, pal$red),
         legend = c("European call payoff",
                    "Replicating portfolio (Delta*S1 + D*(1+r))",
                    "The two possible outcomes (payoffs match)"))
}
save_fig("fig2_replication_payoff", fig_replication, 9, 5.5)


# ---- Figure 3: Q1(d)(i) arbitrage cash flows per call -------------------------
fig_arbitrage <- function() {
  par(mar = c(4.5, 5, 3, 1))
  flows <- c(`Sell call\n(+15)` = Cmkt,
             `Buy Delta shares\n(-50)` = -Delta * S0,
             `Borrow\n(+36.364)` = -D,
             `Riskless profit\nper call` = Cmkt - C0)
  bp <- barplot(flows, col = c(pal$blue, pal$red, pal$green, pal$orange),
                ylab = "Time-0 cash flow per call ($)",
                main = "Q1(d)(i): Time-0 cash flows of the arbitrage (per call sold)",
                ylim = c(min(flows) * 1.25, max(flows) * 1.25), border = NA, cex.names = 0.85)
  abline(h = 0)
  text(bp, flows, labels = sprintf("%.3f", flows),
       pos = ifelse(flows >= 0, 3, 1), font = 2, cex = 0.95)
  mtext(sprintf("Sell %s calls => profit = $%s",
                format(round(n_calls), big.mark = ","),
                format(round(n_calls * (Cmkt - C0)), big.mark = ",")),
        side = 3, line = 0, cex = 0.9)
}
save_fig("fig3_arbitrage_cashflows", fig_arbitrage, 9, 5.5)


# ---- Figure 4: Q2(a)-(c) probabilities and option prices vs u ------------------
fig_q2_prices <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  ug <- seq(1.001, 2, length.out = 300)
  # Left: q and 1 - q
  plot(ug, q_fun(ug), type = "l", lwd = 3, col = pal$blue, ylim = c(0, 1),
       xlab = "u", ylab = "Risk-neutral probability",
       main = "Q2(a): q = 1/(u+1)")
  lines(ug, 1 - q_fun(ug), lwd = 3, col = pal$red)
  abline(h = 0.5, lty = 3, col = pal$grey)
  legend("right", bty = "n", lwd = 3, col = c(pal$blue, pal$red),
         legend = c("q  (up)", "1 - q  (down)"))
  # Right: option prices
  plot(ug, call_CADAUD(ug), type = "l", lwd = 3, col = pal$orange,
       xlab = "u", ylab = "Time-0 price of ATM call",
       main = "Q2(b),(c): ATM call prices (S0 = 1, K = 1)")
  lines(ug, call_AUDCAD(ug), lwd = 3, col = pal$blue)
  legend("topleft", bty = "n", lwd = 3, col = c(pal$blue, pal$orange),
         legend = c("AUD/CAD: (u-1)/(u+1)", "CAD/AUD: u(u-1)/(u+1)"))
}
save_fig("fig4_q2_prices", fig_q2_prices, 10, 5)


# ---- Figure 5: Q2(d) simulated paths of S and Y = 1/S --------------------------
set.seed(3150)
sig_p <- 0.4; n_steps <- 250; n_sim <- 20000
tgrid <- seq(0, 1, length.out = n_steps + 1)
dWm   <- matrix(rnorm(n_sim * n_steps, 0, sqrt(1 / n_steps)), n_sim, n_steps)
Wm    <- cbind(0, t(apply(dWm, 1, cumsum)))
Sm    <- exp(sig_p * Wm - 0.5 * sig_p^2 * matrix(tgrid, n_sim, n_steps + 1, byrow = TRUE))
Ym    <- 1 / Sm

fig_paths <- function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  k <- 40
  matplot(tgrid, t(Sm[1:k, ]), type = "l", lty = 1, lwd = 0.8, col = pal$lightgrey,
          xlab = "t", ylab = "S(t)", main = "AUD/CAD: dS = sigma S dW",
          ylim = range(Sm[1:k, ]))
  lines(tgrid, colMeans(Sm), lwd = 3, col = pal$blue)
  abline(h = 1, lty = 3, col = pal$grey)
  legend("topleft", bty = "n", lwd = c(1, 3), col = c(pal$grey, pal$blue),
         legend = c("Sample paths", "Monte Carlo mean (= 1: martingale)"))
  matplot(tgrid, t(Ym[1:k, ]), type = "l", lty = 1, lwd = 0.8, col = pal$lightgrey,
          xlab = "t", ylab = "Y(t) = 1/S(t)",
          main = expression("CAD/AUD: dY = " * sigma^2 * Y * dt - sigma * Y * dW),
          ylim = range(Ym[1:k, ]))
  lines(tgrid, colMeans(Ym), lwd = 3, col = pal$orange)
  lines(tgrid, exp(sig_p^2 * tgrid), lwd = 2, lty = 2, col = pal$red)
  legend("topleft", bty = "n", lwd = c(1, 3, 2), lty = c(1, 1, 2),
         col = c(pal$grey, pal$orange, pal$red),
         legend = c("Sample paths", "Monte Carlo mean",
                    expression("Theory: " * exp(sigma^2 * t))))
}
save_fig("fig5_q2d_sample_paths", fig_paths, 11, 5.2)


# ---- Figure 6: Q2(e) flash crash and the behaviour of Y ------------------------
set.seed(1)
n_f <- 1000; dtf <- 1 / n_f; tf <- seq(0, 1, length.out = n_f + 1)
sig_f <- 0.2
Sf <- exp(cumsum(c(0, -0.5 * sig_f^2 * dtf + sig_f * rnorm(n_f, 0, sqrt(dtf)))))
i_crash <- which(tf >= 0.6)[1]
Sf[i_crash:(n_f + 1)] <- Sf[i_crash:(n_f + 1)] *
  exp(seq(0, -4, length.out = n_f + 2 - i_crash))      # S falls rapidly towards 0
Yf <- 1 / Sf

fig_flash <- function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.7, 3, 1))
  plot(tf, Sf, type = "l", lwd = 2, col = pal$blue, xlab = "t", ylab = "S(t)",
       main = "AUD/CAD with a flash crash")
  abline(v = tf[i_crash], lty = 3, col = pal$red)
  plot(tf, Yf, type = "l", lwd = 2, col = pal$orange, xlab = "t", ylab = "Y(t) = 1/S(t)",
       main = "CAD/AUD explodes")
  abline(v = tf[i_crash], lty = 3, col = pal$red)
  plot(tf, sig_f^2 * Yf, type = "l", lwd = 2, col = pal$red, log = "y",
       xlab = "t", ylab = "Coefficient size (log scale)",
       main = "Drift and volatility of Y over time",
       ylim = range(c(sig_f^2 * Yf, sig_f * Yf)))
  lines(tf, sig_f * Yf, lwd = 2, col = pal$green)
  legend("topleft", bty = "n", lwd = 2, col = c(pal$red, pal$green),
         legend = c(expression("Drift  " * sigma^2 * Y), expression("Abs. vol  " * sigma * Y)))
  Sg <- exp(seq(log(1), log(0.001), length.out = 200))
  plot(Sg, sig_f^2 / Sg, type = "l", lwd = 2, col = pal$red, log = "xy",
       xlab = "S(t)  (log scale)", ylab = "Coefficient size (log scale)",
       main = "As S -> 0, both blow up",
       ylim = range(c(sig_f^2 / Sg, sig_f / Sg)))
  lines(Sg, sig_f / Sg, lwd = 2, col = pal$green)
  legend("topright", bty = "n", lwd = 2, col = c(pal$red, pal$green),
         legend = c(expression(sigma^2 / S), expression(sigma / S)))
}
save_fig("fig6_q2e_flash_crash", fig_flash, 10, 7.5)

cat("\nDone. All figures are in the ./figures folder.\n")

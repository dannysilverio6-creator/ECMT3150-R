# ECMT3150 Financial Econometrics - Assignment 1
# Danny Silverio, SID 550671578
#
# Q1  TSLA momentum strategies, Jul 2022 - Jun 2026
# Q2  Australian unemployment rate, ARIMA
# Q3  y_t = a + b*y_{t-4} + e_t
#
# Footnotes on the question sheet: volatility = sd of log returns,
# Sharpe = mean log return / volatility with rf = 0. Answers to 3 dp.

rm(list = ls())
graphics.off()

library(quantmod)
library(xts)
library(zoo)
library(forecast)
library(tseries)
library(urca)
library(lmtest)

options(width = 100, scipen = 6)
set.seed(3150)

TD <- 252            # trading days per year
VW <- 20             # vol window for the part (e) rule
SAVEPLOTS <- TRUE
if (SAVEPLOTS) dir.create("fig", showWarnings = FALSE)

f3 <- function(x) formatC(x, format = "f", digits = 3)
pv <- function(p) if (p < 0.001) "< 0.001" else paste("=", f3(p))
hd <- function(s) cat("\n", strrep("=", 70), "\n ", s, "\n", strrep("=", 70), "\n", sep = "")

# data files should sit in the same folder as this script - if not, look around
findfile <- function(nm, pat) {
  if (file.exists(nm)) return(nm)
  d <- c(".", path.expand(c("~", "~/Downloads", "~/Desktop", "~/Documents")))
  hit <- list.files(d[dir.exists(d)], pattern = pat, full.names = TRUE, ignore.case = TRUE)
  if (length(hit)) return(hit[1])
  stop("Can't find ", nm, ". Set the working directory to the folder holding it ",
       "(Session > Set Working Directory > To Source File Location).", call. = FALSE)
}

png2 <- function(f, w = 1900, h = 1200) if (SAVEPLOTS) png(file.path("fig", f), w, h, res = 200)
off  <- function() if (SAVEPLOTS) invisible(dev.off())


###############################################################################
# QUESTION 1
###############################################################################

## ---- (a) data and plots -----------------------------------------------------
hd("Q1 (a)")

px <- tryCatch(getSymbols("TSLA", src = "yahoo", from = "2022-07-01",
                          to = "2026-07-01", auto.assign = FALSE),
               error = function(e) NULL)
if (is.null(px)) {                      # fall back to the saved copy
  raw <- read.csv(findfile("TSLA_2022-07_2026-06.csv", "TSLA.*csv"), stringsAsFactors = FALSE)
  px  <- xts(as.matrix(raw[, -1]), order.by = as.Date(raw$Date))
}
px <- px["2022-07-01/2026-06-30"]

P  <- Cl(px)                            # question asks for the close price
p  <- log(P)
r  <- na.omit(diff(p))
Rs <- exp(r) - 1                        # simple returns, for compounding
rn <- as.numeric(r); Rv <- as.numeric(Rs); n <- length(rn); dt <- index(r)

# close vs adjusted - TSLA pays no dividends so these should be identical
cat("max |Adj/Close - 1| =", max(abs(as.numeric(Ad(px))/as.numeric(Cl(px)) - 1)), "\n")
cat("prices:", NROW(P), " returns:", n, " first/last:",
    round(as.numeric(first(P)), 2), "/", round(as.numeric(last(P)), 2), "\n\n")

sk <- mean((rn - mean(rn))^3)/sd(rn)^3
ku <- mean((rn - mean(rn))^4)/sd(rn)^4
jb <- n/6 * (sk^2 + (ku - 3)^2/4)

cat("mean          ", f3(100*mean(rn)), "% per day  (", f3(100*mean(rn)*TD), "% pa)\n")
cat("sd            ", f3(100*sd(rn)),   "% per day  (", f3(100*sd(rn)*sqrt(TD)), "% pa)\n")
cat("min / max     ", f3(100*min(rn)), "% /", f3(100*max(rn)), "%\n")
cat("skew          ", f3(sk), "\n")
cat("excess kurt   ", f3(ku - 3), "\n")
cat("Jarque-Bera   ", f3(jb), " p", pv(pchisq(jb, 2, lower.tail = FALSE)), "\n")

lb1 <- Box.test(rn, 10, "Ljung-Box"); lb2 <- Box.test(rn^2, 10, "Ljung-Box")
cat("LB(10) r_t    ", f3(lb1$statistic), " p", pv(lb1$p.value), "\n")
cat("LB(10) r_t^2  ", f3(lb2$statistic), " p", pv(lb2$p.value), "\n")

# Engle ARCH LM test - regress squared demeaned returns on 10 of their own lags
e2 <- (rn - mean(rn))^2
X  <- embed(e2, 11)
LM <- nrow(X) * summary(lm(X[, 1] ~ X[, -1]))$r.squared
cat("ARCH LM(10)   ", f3(LM), " p", pv(pchisq(LM, 10, lower.tail = FALSE)), "\n")

png2("Q1a_logprice_logreturns.png", 1900, 1500)
par(mfrow = c(2,1), mar = c(3.2,4.4,2.6,1.2), mgp = c(2.5,0.8,0))
plot(index(p), as.numeric(p), type = "l", col = "navy", lwd = 1.7, xlab = "",
     ylab = expression(p[t] == log(P[t])), main = "TSLA: daily log-price, Jul 2022 - Jun 2026")
grid(col = "grey88")
plot(dt, rn, type = "l", col = "grey35", lwd = 0.7, xlab = "",
     ylab = expression(r[t] == p[t] - p[t-1]), main = "TSLA: daily log-returns")
abline(h = 0, col = "red", lty = 2); grid(col = "grey88")
off(); par(mfrow = c(1,1))

png2("Q1a_acf.png", 1900, 900)
par(mfrow = c(1,2), mar = c(3.6,4,2.8,1), mgp = c(2.4,0.8,0))
acf(rn,   lag.max = 30, main = expression(paste("ACF of ", r[t], ": no mean predictability")))
acf(rn^2, lag.max = 30, main = expression(paste("ACF of ", r[t]^2, ": volatility clustering")))
off(); par(mfrow = c(1,1))


## ---- (b) simple momentum ----------------------------------------------------
hd("Q1 (b)")

# w_t = 1 if r_{t-1} > 0. lag.xts puts r_{t-1} at date t so there's no look-ahead
rlag  <- as.numeric(lag.xts(r, 1))
w_mom <- ifelse(!is.na(rlag) & rlag > 0, 1, 0)
w_bh  <- rep(1, n)

lr <- function(w) log(1 + w*Rv)                       # strategy log returns
sr <- function(w) { x <- lr(w); mean(x)/sd(x) }       # daily Sharpe, rf = 0
wl <- function(w) cumprod(1 + w*Rv)                   # wealth path
dd <- function(W) min(W/cummax(W) - 1)
tn <- function(w) mean(c(abs(w[1]), abs(diff(w))))*TD # annualised turnover

stats <- function(w) {
  x <- lr(w); W <- wl(w)
  c(mean = mean(x), sd = sd(x), mean_a = mean(x)*TD, vol_a = sd(x)*sqrt(TD),
    SR = mean(x)/sd(x), SR_a = mean(x)/sd(x)*sqrt(TD), final = W[n],
    CAGR = W[n]^(TD/n) - 1, maxDD = dd(W), turn = tn(w), expo = mean(w))
}
s_mom <- stats(w_mom); s_bh <- stats(w_bh)

cat("days in market:", sum(w_mom), "of", n, "(", f3(100*mean(w_mom)), "%)\n\n")
print(round(rbind(Momentum = s_mom, BuyHold = s_bh), 6))

cat("\nvol ratio mom/bh:", f3(s_mom["sd"]/s_bh["sd"]), " vs 1/sqrt(2) =", f3(1/sqrt(2)), "\n")

# does the sign of yesterday's return actually predict anything?
up <- w_mom == 1
tt <- t.test(rn[up], rn[!up])
cat("acf(1)          ", f3(acf(rn, 5, plot = FALSE)$acf[2]), "\n")
cat("mean r after up ", f3(100*mean(rn[up])),  "%  (n =", sum(up), ")\n")
cat("mean r after dn ", f3(100*mean(rn[!up])), "%  (n =", sum(!up), ")\n")
cat("difference      t =", f3(tt$statistic), " p =", f3(tt$p.value), "\n")

# any 0/1 rule only helps if cov(w, r) > 0
cat("cov(w, r)       ", round(mean(w_mom*rn) - mean(w_mom)*mean(rn), 6), "\n")


## ---- (c) cumulative return on $1 --------------------------------------------
hd("Q1 (c)")

W_mom <- wl(w_mom); W_bh <- wl(w_bh)
cat("log vs simple compounding, max diff:", max(abs(W_mom - exp(cumsum(w_mom*rn)))), "\n\n")
cat("momentum   $", f3(W_mom[n]), "  total", f3(100*(W_mom[n]-1)), "%  CAGR", f3(100*s_mom["CAGR"]), "%  maxDD", f3(100*s_mom["maxDD"]), "%\n")
cat("buy-hold   $", f3(W_bh[n]),  "  total", f3(100*(W_bh[n]-1)),  "%  CAGR", f3(100*s_bh["CAGR"]),  "%  maxDD", f3(100*s_bh["maxDD"]),  "%\n")

png2("Q1c_cumulative_returns.png")
par(mar = c(3.2,4.4,2.8,1.2), mgp = c(2.6,0.8,0))
plot(dt, W_bh, type = "l", col = "firebrick", lwd = 2, ylim = range(c(W_mom, W_bh)),
     xlab = "", ylab = "Value of $1 invested",
     main = "Cumulative return on $1: momentum vs buy-and-hold (TSLA)")
lines(dt, W_mom, col = "navy", lwd = 2)
abline(h = 1, col = "grey45", lty = 3); grid(col = "grey88")
legend("topleft", bty = "n", lwd = 2, col = c("navy", "firebrick"),
       legend = c(paste0("Momentum     (final $", f3(W_mom[n]), ")"),
                  paste0("Buy-and-hold (final $", f3(W_bh[n]), ")")))
off()


## ---- (d) Sharpe ratios ------------------------------------------------------
hd("Q1 (d)")

cat("momentum      ", f3(s_mom["SR"]), " per day  (ann", f3(s_mom["SR_a"]), ")\n")
cat("buy-and-hold  ", f3(s_bh["SR"]),  " per day  (ann", f3(s_bh["SR_a"]),  ")\n")
cat("difference    ", f3(s_mom["SR"] - s_bh["SR"]), "\n\n")

# Lo (2002) standard error, iid case
se <- function(S) sqrt((1 + 0.5*S^2)/n)
cat("se of daily SR ~", f3(se(s_bh["SR"])), " (ann", f3(se(s_bh["SR"])*sqrt(TD)), ")\n")

# Jobson-Korkie test of equal Sharpe ratios, Memmel (2002) correction
memmel <- function(x, y) {
  mx <- mean(x); my <- mean(y); sx <- sd(x); sy <- sd(y); rho <- cor(x, y)
  v <- (2 - 2*rho + 0.5*(mx/sx)^2 + 0.5*(my/sy)^2 - (mx*my)/(sx*sy)*(1 + rho^2))/length(x)
  z <- (my*sx - mx*sy)/(sx*sy)/sqrt(v)
  c(z = z, p = 2*pnorm(-abs(z)))
}
m <- memmel(lr(w_bh), lr(w_mom))
cat("Memmel test    z =", f3(m["z"]), " p =", f3(m["p"]), "\n")
tp <- t.test(lr(w_mom) - lr(w_bh))
cat("paired t-test  t =", f3(tp$statistic), " p =", f3(tp$p.value), "\n")

# net of transaction costs, charged on turnover
srnet <- function(w, bp) {
  tu <- c(abs(w[1]), abs(diff(w)))
  x <- log(1 + w*Rv - (bp/1e4)*tu); mean(x)/sd(x)
}
cat("\nSharpe net of costs:\n")
ct <- t(sapply(c(0,1,2,5,10,20), function(b) c(BuyHold = srnet(w_bh, b), Momentum = srnet(w_mom, b))))
rownames(ct) <- paste0(c(0,1,2,5,10,20), " bp"); print(round(ct, 4))

be <- uniroot(function(b) srnet(w_mom, b) - s_bh["SR"], c(0, 40))$root
cat("break-even cost:", f3(be), "bp\n")
cat("turnover: momentum", f3(s_mom["turn"]), " buy-hold", f3(s_bh["turn"]), "\n")


## ---- (e) proposed strategy --------------------------------------------------
hd("Q1 (e)")

# w_t = min(1, max(0, r_{t-1}/sigma_{t-1})) - part (b) keeps only the sign of
# this, the proposal keeps the size as well
sig   <- as.numeric(lag.xts(rollapplyr(r, VW, sd, fill = NA), 1))
z     <- rlag/sig; z[!is.finite(z)] <- 0
w_new <- pmin(1, pmax(0, z)); w_new[is.na(w_new)] <- 0
s_new <- stats(w_new); W_new <- wl(w_new)

# check for look-ahead: rebuild every weight from t-1 info only
bad <- 0
for (t in (VW+5):n) {
  wt <- min(1, max(0, rn[t-1]/sd(rn[(t-VW):(t-1)])))
  if (abs(wt - w_new[t]) > 1e-12) bad <- bad + 1
}
cat("look-ahead check:", bad, "of", n-VW-4, "weights disagree\n\n")

print(round(rbind(BuyHold = s_bh, Momentum = s_mom, Proposed = s_new), 6))
cat("\ncov(w,r): momentum", round(mean(w_mom*rn) - mean(w_mom)*mean(rn), 6),
    " proposed", round(mean(w_new*rn) - mean(w_new)*mean(rn), 6), "\n")

# where does the gain come from? middle row is vol targeting on the binary signal
sc     <- pmin(1, (0.30/sqrt(TD))/sig); sc[is.na(sc)] <- 0
w_size <- w_mom*sc
cat("\nattribution:\n")
att <- rbind(BuyHold = stats(w_bh), SignOnly = stats(w_mom),
             SignXVolTarget = stats(w_size), Proposed = stats(w_new))
print(round(att[, c("SR","SR_a","mean_a","vol_a","final","maxDD","turn","expo")], 4))

# conditional means by quintile of the standardised signal
ok <- which(is.finite(z) & !is.na(sig))
q  <- cut(z[ok], quantile(z[ok], 0:5/5), include.lowest = TRUE, labels = paste0("Q", 1:5))
cat("\nmean r_t by quintile of z_t (annualised %):\n")
print(round(t(sapply(levels(q), function(g) {
  i <- ok[q == g]
  c(n = length(i), zmin = min(z[i]), zmax = max(z[i]),
    ann.mean = 100*mean(rn[i])*TD, ann.sd = 100*sd(rn[i])*sqrt(TD),
    se = 100*TD*sd(rn[i])/sqrt(length(i)))
})), 3))
pr <- summary(lm(rn[ok] ~ pmax(0, z[ok])))
cat("regression r_t on max(0,z_t): b =", f3(coef(pr)[2,1]),
    " t =", f3(coef(pr)[2,3]), " p =", f3(coef(pr)[2,4]), "\n")

# robustness 1: does it depend on the parameters?
mk <- function(lam, vw, cap) {
  s <- as.numeric(lag.xts(rollapplyr(r, vw, sd, fill = NA), 1))
  zz <- rlag/(lam*s); zz[!is.finite(zz)] <- 0
  ww <- pmin(cap, pmax(0, zz)); ww[is.na(ww)] <- 0; ww
}
g <- expand.grid(lam = c(0.5,0.75,1,1.5,2,3), vw = c(10,20,40,60,120), cap = c(0.5,1,1.5,2))
mid <- floor(n/2); i1 <- 1:mid; i2 <- (mid+1):n
sub <- function(w, i) { x <- log(1 + w[i]*Rv[i]); mean(x)/sd(x) }
for (k in seq_len(nrow(g))) {
  ww <- mk(g$lam[k], g$vw[k], g$cap[k])
  g$SR[k] <- sr(ww); g$SR5[k] <- srnet(ww, 5)
  g$train[k] <- sub(ww, i1); g$test[k] <- sub(ww, i2)
}
bar <- max(s_mom["SR"], s_bh["SR"])
cat("\ngrid of", nrow(g), "parameter combos, bar to beat =", f3(bar), "\n")
cat("  beat both       :", f3(100*mean(g$SR > bar)), "%\n")
cat("  beat both @ 5bp :", f3(100*mean(g$SR5 > bar)), "%\n")
cat("  test half beats :", f3(100*mean(g$test > bar)), "%\n")
cat("  min/median/max  :", f3(min(g$SR)), "/", f3(median(g$SR)), "/", f3(max(g$SR)), "\n")
cat("  baseline is at the", f3(100*mean(g$SR <= s_new["SR"])), "th percentile\n")

# robustness 2: split sample
cat("\nsplit sample (train", format(dt[1]), "-", format(dt[mid]),
    "| test", format(dt[mid+1]), "-", format(dt[n]), "):\n")
oos <- rbind(Train = c(sub(w_bh,i1), sub(w_mom,i1), sub(w_new,i1)),
             Test  = c(sub(w_bh,i2), sub(w_mom,i2), sub(w_new,i2)))
colnames(oos) <- c("BuyHold","Momentum","Proposed"); print(round(oos, 4))

# robustness 3: costs, and dropping the best days
cat("\nSharpe net of costs:\n")
ct2 <- t(sapply(c(0,1,2,5,10,20), function(b)
  sapply(list(w_bh, w_mom, w_new), srnet, bp = b)))
dimnames(ct2) <- list(paste0(c(0,1,2,5,10,20)," bp"), c("BuyHold","Momentum","Proposed"))
print(round(ct2, 4))

cat("\ndropping the proposed strategy's best days:\n")
for (k in c(0,5,10,20)) {
  keep <- setdiff(1:n, order(lr(w_new), decreasing = TRUE)[seq_len(k)])
  f <- function(w) { x <- lr(w)[keep]; f3(mean(x)/sd(x)) }
  cat("  drop", sprintf("%2d", k), ":", f(w_bh), f(w_mom), f(w_new), "\n")
}

cat("\nMemmel: proposed vs momentum z =", f3(memmel(lr(w_mom), lr(w_new))["z"]),
    " p =", f3(memmel(lr(w_mom), lr(w_new))["p"]),
    "| vs buy-hold z =", f3(memmel(lr(w_bh), lr(w_new))["z"]),
    " p =", f3(memmel(lr(w_bh), lr(w_new))["p"]), "\n")

png2("Q1e_strategy_comparison.png")
par(mar = c(3.2,4.4,2.8,1.2), mgp = c(2.6,0.8,0))
plot(dt, W_bh, type = "l", col = "firebrick", lwd = 2, ylim = range(c(W_bh,W_mom,W_new)),
     xlab = "", ylab = "Value of $1 invested", main = "Part (e): cumulative return on $1")
lines(dt, W_mom, col = "navy", lwd = 2); lines(dt, W_new, col = "darkgreen", lwd = 2.2)
abline(h = 1, col = "grey45", lty = 3); grid(col = "grey88")
legend("topleft", bty = "n", lwd = 2, cex = 0.9,
       col = c("darkgreen","navy","firebrick"),
       legend = c(paste("Proposed     SR", f3(s_new["SR"]), "  final $", f3(W_new[n])),
                  paste("Momentum     SR", f3(s_mom["SR"]), "  final $", f3(W_mom[n])),
                  paste("Buy-and-hold SR", f3(s_bh["SR"]),  "  final $", f3(W_bh[n]))))
off()

png2("Q1e_weights_and_vol.png", 1900, 1500)
par(mfrow = c(2,1), mar = c(3.2,4.4,2.6,1.2), mgp = c(2.5,0.8,0))
plot(dt, sig*sqrt(TD), type = "l", col = "grey35", lwd = 0.9, xlab = "",
     ylab = "Annualised volatility", main = "Rolling 20-day volatility of TSLA, lagged one day")
grid(col = "grey88")
plot(dt, w_new, type = "h", col = "darkgreen", lwd = 0.8, xlab = "", ylab = "Weight in TSLA",
     main = "Position taken by the proposed strategy (zero on half of all days)")
abline(h = c(0,1), col = "grey55", lty = 3); grid(col = "grey88")
off(); par(mfrow = c(1,1))

png2("Q1e_robustness.png", 1900, 900)
par(mfrow = c(1,2), mar = c(3.8,4,2.8,1), mgp = c(2.4,0.8,0))
hist(g$SR, breaks = 20, col = "grey85", border = "white",
     xlim = range(c(bar, g$SR))*c(0.95, 1.02), xlab = "Daily Sharpe ratio",
     main = paste("All", nrow(g), "parameter combinations"))
abline(v = bar, col = "firebrick", lwd = 2, lty = 2)
abline(v = s_new["SR"], col = "darkgreen", lwd = 2)
legend("topleft", bty = "n", lwd = 2, lty = c(2,1), cex = 0.8,
       col = c("firebrick","darkgreen"),
       legend = c(paste0("bar to beat (", f3(bar), ")"),
                  paste0("baseline (", f3(s_new["SR"]), ")")))
plot(g$train, g$test, pch = 16, col = "darkgreen", cex = 0.7,
     xlab = "Sharpe, training half", ylab = "Sharpe, test half", main = "Split-sample stability")
abline(h = bar, v = bar, col = "firebrick", lty = 2); grid(col = "grey88")
off(); par(mfrow = c(1,1))


###############################################################################
# QUESTION 2
###############################################################################

u  <- read.csv(findfile("aus_unemploy_rate.csv", "unemploy.*csv"),
               fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
y  <- ts(u[[2]], start = c(2006, 5), frequency = 12)
tr <- window(y, end = c(2019, 12))

png2("Q2a_series.png", 1900, 1000)
par(mar = c(3.2,4.2,2.8,1.2), mgp = c(2.5,0.8,0))
plot(y, col = "grey55", lwd = 1.6, xlab = "", ylab = "Unemployment rate (%)",
     main = "Australian monthly unemployment rate, May 2006 - May 2026")
lines(tr, col = "navy", lwd = 2)
abline(v = 2020, col = "firebrick", lty = 2); grid(col = "grey88")
legend("topleft", bty = "n", lwd = 2, lty = c(1,1,2), col = c("navy","grey55","firebrick"),
       legend = c("Estimation sample (to Dec 2019)", "Held out (Jan 2020 on)", "January 2020"))
off()


## ---- (a) ADF test -----------------------------------------------------------
hd("Q2 (a)")

# H0: gamma = 0 (unit root)  vs  H1: gamma < 0 (stationary), 5% level.
# Intercept but no trend - the series has no sustained trend over this window.
a1 <- ur.df(tr, type = "drift", selectlags = "AIC")
print(round(summary(a1)@testreg$coefficients, 6))

a2 <- ur.df(tr, type = "trend", selectlags = "AIC")
a3 <- ur.df(diff(tr), type = "drift", selectlags = "AIC")
adf <- rbind(`level, drift` = c(a1@teststat[1], a1@cval[1,]),
             `level, trend` = c(a2@teststat[1], a2@cval[1,]),
             `diff,  drift` = c(a3@teststat[1], a3@cval[1,]))
colnames(adf)[1] <- "stat"
cat("\n"); print(round(adf, 3))
cat("\nlevel: ", f3(a1@teststat[1]), " > ", f3(a1@cval[1,"5pct"]), " -> cannot reject unit root\n", sep = "")
cat("diff : ", f3(a3@teststat[1]), " < ", f3(a3@cval[1,"5pct"]), " -> reject, so u_t is I(1), d = 1\n", sep = "")

print(adf.test(tr))
k1 <- suppressWarnings(kpss.test(tr)); k2 <- suppressWarnings(kpss.test(diff(tr)))
cat("KPSS level:", f3(k1$statistic), " p", pv(k1$p.value),
    " | diff:", f3(k2$statistic), " p =", f3(k2$p.value), "\n")
cat("ndiffs:", ndiffs(tr, test = "adf"), ndiffs(tr, test = "kpss"), " nsdiffs:", nsdiffs(tr), "\n")

png2("Q2a_acf.png", 1900, 1560)
par(mfrow = c(2,2), mar = c(3.6,4,3.4,1), mgp = c(2.3,0.8,0))
acf (tr,       36, main = "ACF: level");            pacf(tr,       36, main = "PACF: level")
acf (diff(tr), 36, main = "ACF: first difference"); pacf(diff(tr), 36, main = "PACF: first difference")
off(); par(mfrow = c(1,1))


## ---- (b) fit an ARIMA -------------------------------------------------------
hd("Q2 (b)")

roots <- function(f) {
  ar <- if (length(f$model$phi)   && any(f$model$phi   != 0)) 1/Mod(polyroot(c(1, -f$model$phi)))   else numeric(0)
  ma <- if (length(f$model$theta) && any(f$model$theta != 0)) 1/Mod(polyroot(c(1,  f$model$theta))) else numeric(0)
  list(ar = ar, ma = ma, ok = (!length(ma) || all(ma < 1 - 1e-8)) && (!length(ar) || all(ar < 1 - 1e-8)))
}

gr <- data.frame()
for (pp in 0:4) for (qq in 0:4) {
  f <- tryCatch(Arima(tr, order = c(pp,1,qq), method = "ML"), error = function(e) NULL)
  if (!is.null(f)) gr <- rbind(gr, data.frame(p = pp, q = qq, AIC = AIC(f), AICc = f$aicc,
    BIC = BIC(f), sigma2 = f$sigma2, ok = roots(f)$ok,
    LB = Box.test(residuals(f), 24, "Ljung-Box", fitdf = pp+qq)$p.value))
}
cat("top 10 by AIC (d = 1):\n"); print(head(round(gr[order(gr$AIC), ], 4), 10), row.names = FALSE)

# AIC minimum is (2,1,4) but it's not invertible and fails Ljung-Box, so it's out
chk <- Arima(tr, order = c(2,1,4), method = "ML")
cat("\nARIMA(2,1,4): inverse MA roots", f3(roots(chk)$ma),
    "| LB p =", f3(Box.test(residuals(chk), 24, "Ljung-Box", fitdf = 6)$p.value), "\n")

good <- subset(gr, ok & LB > 0.05)
sel  <- good[which.min(good$AIC), ]
fit  <- Arima(tr, order = c(sel$p, 1, sel$q), method = "ML")
cat("\nselected: ARIMA(", sel$p, ",1,", sel$q, ")\n", sep = "")
print(fit); cat("\n"); print(coeftest(fit))
rr <- roots(fit)
cat("\ninverse AR roots", f3(rr$ar), " inverse MA roots", f3(rr$ma), "\n")
lbf <- Box.test(residuals(fit), 24, "Ljung-Box", fitdf = sel$p + sel$q)
cat("LB(24) p =", f3(lbf$p.value),
    " | LB(24) on squares p", pv(Box.test(residuals(fit)^2, 24, "Ljung-Box")$p.value),
    " | JB p =", f3(jarque.bera.test(residuals(fit))$p.value), "\n")
zr <- polyroot(c(1, -fit$model$phi))
cat("AR roots", round(zr, 3), " -> cycle of", f3(2*pi/abs(Arg(1/zr)[1])), "months\n")

png2("Q2b_residuals.png", 1900, 1560)
par(mfrow = c(2,2), mar = c(3.6,4,3.2,1), mgp = c(2.3,0.8,0))
plot(residuals(fit), col = "navy", xlab = "", ylab = "residual", main = "Residuals over time")
abline(h = 0, col = "grey45"); grid(col = "grey90")
acf(residuals(fit), 24, main = "ACF of the residuals")
hist(residuals(fit), breaks = 25, col = "grey85", border = "white", freq = FALSE,
     xlab = "residual", main = "Residual distribution")
curve(dnorm(x, 0, sd(residuals(fit))), add = TRUE, col = "firebrick", lwd = 2)
qqnorm(residuals(fit), pch = 16, cex = 0.5, col = "grey35"); qqline(residuals(fit), col = "firebrick", lwd = 2)
off(); par(mfrow = c(1,1))


## ---- (c) forecast Jan-Jun 2020 ----------------------------------------------
hd("Q2 (c)")

fc  <- forecast(fit, h = 6, level = 95)
act <- window(y, start = c(2020,1), end = c(2020,6))
out <- data.frame(month = c("Jan","Feb","Mar","Apr","May","Jun"),
                  fcast = as.numeric(fc$mean),
                  se    = as.numeric((fc$upper - fc$mean)/qnorm(0.975)),
                  lo    = as.numeric(fc$lower), hi = as.numeric(fc$upper),
                  actual = as.numeric(act))
out$z      <- (out$actual - out$fcast)/out$se
out$inside <- out$actual >= out$lo & out$actual <= out$hi
print(cbind(out[1], round(out[-1][-6], 3), inside = out$inside), row.names = FALSE)
cat("\nRMSE over the 6 months:", f3(sqrt(mean((out$actual - out$fcast)^2))), "\n")

png2("Q2c_forecast.png", 1900, 1100)
par(mar = c(3.2,4.2,2.8,1.2), mgp = c(2.5,0.8,0))
plot(fc, xlim = c(2015, 2020.6), ylim = c(3.8, 7.8), xlab = "", ylab = "Unemployment rate (%)",
     main = paste0("ARIMA(", sel$p, ",1,", sel$q, ") forecast with 95% interval"),
     fcol = "navy", shadecols = "#C6D2EC", flwd = 2)
lines(window(y, start = c(2015,1), end = c(2020,6)), col = "grey40", lwd = 1.4)
points(act, col = "firebrick", pch = 19, cex = 0.8); lines(act, col = "firebrick", lwd = 2)
grid(col = "grey88")
legend("topleft", bty = "n", lwd = 2, pch = c(NA,19), col = c("navy","firebrick"),
       legend = c("Forecast + 95% interval", "Realised, Jan-Jun 2020"))
off()


## ---- (d) June 2020 interval -------------------------------------------------
hd("Q2 (d)")

cat("June 2020: forecast", f3(out$fcast[6]), " se", f3(out$se[6]),
    " -> 95% CI [", f3(out$lo[6]), ",", f3(out$hi[6]), "]\n")
cat("realised", f3(out$actual[6]), "-> OUTSIDE, ", f3(out$z[6]), "standard errors above\n")
cat("prob of an outcome this extreme:", signif(2*pnorm(-abs(out$z[6])), 3), "\n")
cat("Jan-Mar are inside the interval, Apr-Jun are outside - the break is April 2020.\n")


## ---- (e) refit on the whole sample ------------------------------------------
hd("Q2 (e)")

af <- ur.df(y, type = "drift", selectlags = "AIC")
cat("full-sample ADF:", f3(af@teststat[1]), "vs", f3(af@cval[1,"5pct"]), "-> still d = 1\n\n")

gf <- data.frame()
for (pp in 0:4) for (qq in 0:4) {
  f <- tryCatch(Arima(y, order = c(pp,1,qq), method = "ML"), error = function(e) NULL)
  if (!is.null(f)) gf <- rbind(gf, data.frame(p = pp, q = qq, AIC = AIC(f),
                                              AICc = f$aicc, BIC = BIC(f), sigma2 = f$sigma2))
}
cat("top 8 by AIC:\n"); print(head(round(gf[order(gf$AIC), ], 4), 8), row.names = FALSE)
b2 <- gf[which.min(gf$AIC), ]
cat("\nAIC min:", paste0("ARIMA(", b2$p, ",1,", b2$q, ")"),
    " BIC min:", paste0("ARIMA(", gf$p[which.min(gf$BIC)], ",1,", gf$q[which.min(gf$BIC)], ")"), "\n")

fit_f  <- Arima(y, order = c(b2$p, 1, b2$q), method = "ML")
fit_f2 <- Arima(y, order = c(sel$p, 1, sel$q), method = "ML")   # part (b) order, forced
print(fit_f)
cat("\npart (b) order forced onto the full sample:\n"); print(coeftest(fit_f2))
cat("AIC", f3(AIC(fit_f2)), "vs random walk", f3(AIC(fit_f)), "\n")

diag <- function(f, k) {
  r <- residuals(f)
  c(sigma2 = f$sigma2, sd = sd(r),
    LB.p = Box.test(r, 24, "Ljung-Box", fitdf = k)$p.value,
    JB = jarque.bera.test(r)$statistic,
    skew = mean((r-mean(r))^3)/sd(r)^3, exkurt = mean((r-mean(r))^4)/sd(r)^4 - 3,
    maxabs = max(abs(r)), n3sd = sum(abs(r) > 3*sd(r)))
}
cat("\nresidual comparison:\n")
print(round(t(rbind(`train (2,1,2)` = diag(fit, sel$p+sel$q),
                    `full (0,1,0)`  = diag(fit_f, 0),
                    `full (2,1,2)`  = diag(fit_f2, sel$p+sel$q))), 4))
cat("\nsigma^2 rises", f3(100*(fit_f$sigma2/fit$sigma2 - 1)), "%\n")

rf <- residuals(fit_f)
cat("largest full-sample residuals:\n")
for (i in order(abs(rf), decreasing = TRUE)[1:4])
  cat("  ", format(time(rf)[i]), ":", f3(rf[i]), "(", f3(rf[i]/sd(rf)), "sd )\n")

# what happens if the break is modelled instead of ignored
covid <- as.numeric(time(y) >= 2020.2 & time(y) <= 2021.0)
fd <- Arima(y, order = c(sel$p,1,sel$q), xreg = covid, method = "ML")
cat("\nwith a COVID dummy:\n"); print(coeftest(fd))
cat("sigma^2", f3(fd$sigma2), " AIC", f3(AIC(Arima(y, order=c(0,1,0), xreg=covid, method="ML"))), "\n")

png2("Q2e_residuals.png", 1900, 1560)
par(mfrow = c(2,2), mar = c(3.6,4,3.4,1), mgp = c(2.3,0.8,0))
plot(residuals(fit), col = "navy", xlab = "", ylab = "residual",
     main = paste0("Residuals: training ARIMA(", sel$p, ",1,", sel$q, ")"))
abline(h = 0, col = "grey45"); grid(col = "grey90")
acf(residuals(fit), 24, main = "ACF, training residuals")
plot(rf, col = "firebrick", xlab = "", ylab = "residual",
     main = paste0("Residuals: full-sample ARIMA(", b2$p, ",1,", b2$q, ")"))
abline(h = 0, col = "grey45"); grid(col = "grey90")
acf(rf, 24, main = "ACF, full-sample residuals")
off(); par(mfrow = c(1,1))


###############################################################################
# QUESTION 3   y_t = a + b*y_{t-4} + e_t
###############################################################################
hd("Q3")

a <- 0.02; b <- 0.4; s2 <- 0.01; sg <- sqrt(s2)
yl <- c(-0.015, -0.010, -0.005, 0.000)      # y77 ... y80

mu <- a/(1 - b)                              # (b) unconditional mean
g0 <- s2/(1 - b^2)                           # (c) unconditional variance
rho  <- sapply(1:8, function(k) if (k %% 4 == 0) b^(k/4) else 0)   # (d) ACF
pacf8 <- c(0, 0, 0, b, 0, 0, 0, 0)                                 # (e) PACF

cat("(a) stationary iff |b| < 1; here |b| =", f3(abs(b)), ", roots have modulus", f3(abs(b)^(-1/4)), "\n")
cat("(b) mu   =", f3(mu), "\n")
cat("(c) var  =", f3(g0), " (sd", f3(sqrt(g0)), ")\n")
cat("(d) ACF  :", f3(rho), "\n")
cat("(e) PACF :", f3(pacf8), "\n")

# (f) forecasts - for l = 1..4 the lag is observed, at l = 5 it's a forecast
fc5 <- numeric(5)
for (l in 1:5) fc5[l] <- a + b*(if (l <= 4) yl[l] else fc5[l-4])
cat("(f) forecasts l=1..5:", f3(fc5), "\n")

# (g) psi weights are b^(j/4) at multiples of 4 and zero elsewhere
psi <- sapply(0:40, function(j) if (j %% 4 == 0) b^(j/4) else 0)
sef <- sapply(1:5, function(l) sqrt(s2*sum(psi[1:l]^2)))
cat("(g) se l=1..5       :", f3(sef), "\n")
cat("(h) limit of forecast:", f3(mu), "\n")
cat("(i) limit of se      :", f3(sqrt(g0)), "\n")

# check the algebra by simulating the process
N <- 1e6; burn <- 1000
e <- rnorm(N + burn, 0, sg); ys <- numeric(N + burn)
for (t in 5:(N + burn)) ys[t] <- a + b*ys[t-4] + e[t]
ys <- ys[-(1:burn)]
cat("\nsimulation check (", format(N, big.mark = ","), "draws):\n")
cat("  mean", f3(mean(ys)), "vs", f3(mu), " | var", f3(var(ys)), "vs", f3(g0), "\n")
cat("  ACF ", f3(acf(ys, 8, plot = FALSE)$acf[2:9]), "\n")
cat("  PACF", f3(pacf(ys, 8, plot = FALSE)$acf[1:8]), "\n")

M <- 2e5; err <- matrix(NA_real_, M, 5)
for (m in 1:M) {
  path <- c(yl, numeric(5))
  for (l in 1:5) path[4+l] <- a + b*path[l] + rnorm(1, 0, sg)
  err[m, ] <- path[5:9] - fc5
}
cat("  simulated se", f3(apply(err, 2, sd)), "vs theory", f3(sef), "\n")

png2("Q3_acf_pacf.png", 1900, 900)
par(mfrow = c(1,2), mar = c(3.6,4,2.8,1), mgp = c(2.3,0.8,0))
barplot(rho, names.arg = 1:8, col = "navy", border = NA, ylim = c(-0.05, 0.5),
        xlab = "lag", ylab = expression(rho[k]), main = "Theoretical ACF")
abline(h = 0, col = "grey40"); box()
barplot(pacf8, names.arg = 1:8, col = "darkgreen", border = NA, ylim = c(-0.05, 0.5),
        xlab = "lag", ylab = expression(phi[kk]), main = "Theoretical PACF")
abline(h = 0, col = "grey40"); box()
off(); par(mfrow = c(1,1))

# forecast converging to the unconditional distribution
L <- 40; fl <- numeric(L); hh <- c(yl, numeric(L))
for (l in 1:L) { hh[4+l] <- a + b*hh[l]; fl[l] <- hh[4+l] }
sl <- sapply(1:L, function(l) sqrt(s2*sum(psi[1:min(l, length(psi))]^2)))
png2("Q3_forecast_convergence.png", 1900, 950)
par(mar = c(3.6,4.2,2.8,1.2), mgp = c(2.4,0.8,0))
plot(1:L, fl, type = "o", pch = 16, cex = 0.6, col = "navy", lwd = 2,
     ylim = range(c(fl - 2*sl, fl + 2*sl)), xlab = "forecast horizon  l",
     ylab = "point forecast  y_80(l)",
     main = "Q3: point forecast and 95% band converge to the unconditional distribution")
polygon(c(1:L, rev(1:L)), c(fl + 1.96*sl, rev(fl - 1.96*sl)), col = "#C6D2EC80", border = NA)
lines(1:L, fl, col = "navy", lwd = 2)
abline(h = mu, col = "firebrick", lwd = 2, lty = 2)
abline(h = mu + c(-1.96, 1.96)*sqrt(g0), col = "darkgreen", lty = 3, lwd = 1.6)
grid(col = "grey88")
legend("bottomright", bg = "white", box.col = "grey70", lwd = 2, lty = c(1,2,3), cex = 0.85,
       col = c("navy","firebrick","darkgreen"),
       legend = c("point forecast", paste("unconditional mean", f3(mu)),
                  paste("+/- 1.96 sqrt(g0) =", f3(1.96*sqrt(g0)))))
off()


###############################################################################
# summary of the answers
###############################################################################
hd("Answers")

cat("Q1\n")
q1 <- rbind(`mean daily %` = 100*c(s_bh["mean"], s_mom["mean"], s_new["mean"]),
                  `daily sd %`   = 100*c(s_bh["sd"], s_mom["sd"], s_new["sd"]),
                  `SHARPE (day)` = c(s_bh["SR"], s_mom["SR"], s_new["SR"]),
                  `Sharpe (ann)` = c(s_bh["SR_a"], s_mom["SR_a"], s_new["SR_a"]),
                  `final $1`     = c(s_bh["final"], s_mom["final"], s_new["final"]),
                  `max DD %`     = 100*c(s_bh["maxDD"], s_mom["maxDD"], s_new["maxDD"]),
            `turnover`     = c(s_bh["turn"], s_mom["turn"], s_new["turn"]))
colnames(q1) <- c("BuyHold", "Momentum", "Proposed")
print(round(q1, 3))

cat("\nQ2\n")
cat(" (a) ADF", f3(a1@teststat[1]), "vs", f3(a1@cval[1,"5pct"]), "-> I(1), d = 1\n")
cat(" (b) ARIMA(", sel$p, ",1,", sel$q, "), AIC ", f3(AIC(fit)), ", sigma^2 ", f3(fit$sigma2), "\n", sep = "")
cat(" (c)", f3(out$fcast), "\n")
cat(" (d) [", f3(out$lo[6]), ",", f3(out$hi[6]), "], realised", f3(out$actual[6]), "-> outside\n")
cat(" (e) ARIMA(", b2$p, ",1,", b2$q, "), sigma^2 ", f3(fit$sigma2), " -> ", f3(fit_f$sigma2), "\n", sep = "")

cat("\nQ3\n")
cat(" (a) |b| < 1   (b)", f3(mu), "  (c)", f3(g0), "\n")
cat(" (d)", f3(rho), "\n")
cat(" (e)", f3(pacf8), "\n")
cat(" (f)", f3(fc5), "\n")
cat(" (g)", f3(sef), "\n")
cat(" (h)", f3(mu), "  (i)", f3(sqrt(g0)), "\n")

cat("\n"); print(R.version.string)

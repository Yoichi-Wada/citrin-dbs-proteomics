ci <- function(x, n) { b <- binom.test(x, n)$conf.int; sprintf("%d/%d = %.3f (%.3f-%.3f)", x, n, x/n, b[1], b[2]) }
cat("stage1 sens LoF/LoF:", ci(16,16), "\nstage1 spec (44 non-LoF/LoF):", ci(44,44), "\nstage2 sens LoF/Mis:", ci(6,7), "\noverall sens (23):", ci(22,23), "\noverall spec (37):", ci(26,37), "\nspec DD:", ci(15,21), "\nspec WS:", ci(9,10), "\nspec WW:", ci(2,6), "\nBA overall two-stage: (22/23 + 26/37)/2 =", round((22/23+26/37)/2,3), "\n")
cat("diverse30: overall sens", ci(23,23), " spec", ci(32,37), "\nconsensus30 spec>=0.90 policy: LoF/Mis sens", ci(3,7), " spec", ci(35,37), " overall sens", ci(19,23), "\n")
ppv <- function(se, sp, p) se*p/(se*p + (1-sp)*(1-p)); npv <- function(se, sp, p) sp*(1-p)/(sp*(1-p) + (1-se)*p)
cat("\nPPV/NPV at stated prevalences (two-stage consensus30: se 22/23, sp 26/37; stage 1 alone for LoF/LoF: se 16/16, sp 44/44 with lower CI 0.920):\n")
for (p in c(1/17000, 1/100, 1/10, 23/60)) {
  cat(sprintf("prev %-10s | two-stage PPV %.4f NPV %.5f FP per TP %6.0f | stage1-only (LoF/LoF) PPV %.4f (spec lower CI 0.920 -> PPV %.4f)\n",
              format(p, digits = 3), ppv(22/23, 26/37, p), npv(22/23, 26/37, p), (1 - 26/37) * (1 - p) / ((22/23) * p),
              ppv(1, 1, p), ppv(1, 0.920, p)))
}
cat("\nhypothetical 100% sens / 95% spec at 1/17,000: PPV", sprintf("%.4f", ppv(1, 0.95, 1/17000)),
    " FP per TP", round(0.05 * (1 - 1/17000) / (1/17000)), "\n")

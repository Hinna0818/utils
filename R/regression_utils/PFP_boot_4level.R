#' Calculate Preventable Fraction for the Population (PFP) for 4-level categorical exposure
#'
#' @param model A coxph fitted model including a 4-level factor exposure variable, with level "1" as reference
#' @param exposure_var Character. Name of the 4-level exposure variable in the model
#' @param data The original data.frame used in the model
#' @param n_boot Number of bootstrap resamples (default 1000)
#' @param conf_level Confidence level for CI (default 0.95)
#'
#' @return A list with point estimate, CI lower, upper, and formatted string
#' 
calculate_pfp_boot_4level <- function(model, exposure_var, data, n_boot = 1000, conf_level = 0.95) {
  
  data[[exposure_var]] <- factor(data[[exposure_var]], levels = c("1", "2", "3", "4"))
  
  # 提取 HRs
  coef_table <- summary(model)$coefficients
  name2 <- grep(paste0(exposure_var, "2"), rownames(coef_table), value = TRUE)
  name3 <- grep(paste0(exposure_var, "3"), rownames(coef_table), value = TRUE)
  name4 <- grep(paste0(exposure_var, "4"), rownames(coef_table), value = TRUE)
  
  HR2 <- exp(coef_table[name2, "coef"])
  HR3 <- exp(coef_table[name3, "coef"])
  HR4 <- exp(coef_table[name4, "coef"])
  
  # 主数据下暴露比例
  P2 <- mean(data[[exposure_var]] == "2")
  P3 <- mean(data[[exposure_var]] == "3")
  P4 <- mean(data[[exposure_var]] == "4")
  
  # PFP计算函数
  compute_pfp <- function(p2, p3, p4, hr2, hr3, hr4) {
    num <- p2 * (hr2 - 1) + p3 * (hr3 - 1) + p4 * (hr4 - 1)
    denom <- num + 1
    return(num / denom * 100)
  }
  
  # 主估计值
  pfp_point <- compute_pfp(P2, P3, P4, HR2, HR3, HR4)
  
  # Bootstrap置信区间
  set.seed(2025)
  pfp_boot <- replicate(n_boot, {
    sample_data <- data[sample(nrow(data), replace = TRUE), ]
    p2_b <- mean(sample_data[[exposure_var]] == "2")
    p3_b <- mean(sample_data[[exposure_var]] == "3")
    p4_b <- mean(sample_data[[exposure_var]] == "4")
    compute_pfp(p2_b, p3_b, p4_b, HR2, HR3, HR4)
  })
  
  alpha <- (1 - conf_level) / 2
  ci_lower <- quantile(pfp_boot, probs = alpha)
  ci_upper <- quantile(pfp_boot, probs = 1 - alpha)
  
  return(list(
    estimate = round(pfp_point, 2),
    lower = round(ci_lower, 2),
    upper = round(ci_upper, 2),
    formatted = paste0(round(pfp_point, 2), "% (95% CI: ", round(ci_lower, 2), "% – ", round(ci_upper, 2), "%)")
  ))
}

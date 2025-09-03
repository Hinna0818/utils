#' Calculate Preventable Fraction for the Population (PFP) with Bootstrap Confidence Interval
#'
#' This function estimates the Preventable Fraction for the Population (PFP) based on a fitted Cox proportional hazards model,
#' assuming a 3-level categorical exposure variable (e.g., PsRS = "1", "2", "3") where level "1" is the reference (ideal) group.
#' The function computes the point estimate using hazard ratios (HRs) from the original model and uses bootstrap resampling to
#' estimate a 95% confidence interval by resampling the exposure distribution while keeping the HRs fixed.
#'
#' @param model A fitted \code{coxph} object from the \code{survival} package. The model should include a 3-level factor exposure variable,
#'              with level "1" as the reference group.
#' @param exposure_var Character string. The name of the exposure variable (e.g., "PsRS").
#' @param data A \code{data.frame} containing the variables used in the model, including the exposure variable.
#' @param n_boot Integer. Number of bootstrap replicates to perform. Default is 1000.
#' @param conf_level Numeric. Confidence level for the interval estimate. Default is 0.95.
#'
#'
calculate_pfp_boot <- function(model, exposure_var, data, n_boot = 1000, conf_level = 0.95) {
  
  data[[exposure_var]] <- factor(data[[exposure_var]], levels = c("1", "2", "3"))
  
  # 提取主模型 HR
  coef_table <- summary(model)$coefficients
  name2 <- grep(paste0(exposure_var, "2"), rownames(coef_table), value = TRUE)
  name3 <- grep(paste0(exposure_var, "3"), rownames(coef_table), value = TRUE)
  HR2 <- exp(coef_table[name2, "coef"])
  HR3 <- exp(coef_table[name3, "coef"])
  
  # 主模型下的暴露比例
  P2 <- mean(data[[exposure_var]] == "2")
  P3 <- mean(data[[exposure_var]] == "3")
  
  # PFP公式
  compute_pfp <- function(p2, p3, hr2, hr3) {
    num <- p2 * (hr2 - 1) + p3 * (hr3 - 1)
    denom <- num + 1
    return(num / denom * 100)
  }
  
  # 主估计值
  pfp_point <- compute_pfp(P2, P3, HR2, HR3)
  
  # bootstrap 置信区间
  set.seed(2025) 
  pfp_boot <- replicate(n_boot, {
    sample_data <- data[sample(nrow(data), replace = TRUE), ]
    p2_b <- mean(sample_data[[exposure_var]] == "2")
    p3_b <- mean(sample_data[[exposure_var]] == "3")
    compute_pfp(p2_b, p3_b, HR2, HR3)
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

#' Diabetes Dataset
#'
#' @description
#' A diabetes dataset on 1,000 patients.
#'
#' @docType data
#' @keywords datasets internal
#' @format
#' A data frame with 1,000 rows and 14 variables:
#' \describe{
#'   \item{\code{ID}}{Unique identifier for each patient (unitless).}
#'   \item{\code{No_Pation}}{Patient number (unitless).}
#'   \item{\code{Gender}}{Categorical variable (Female, Male).}
#'   \item{\code{AGE}}{Years (age of the person).}
#'   \item{\code{Urea}}{Chief nitrogenous end product of the metabolic breakdown of proteins in milligrams per deciliter (mg/dL).}
#'   \item{\code{Cr}}{Creatinine ratio (Cr) (mg/dL).}
#'   \item{\code{HbA1c}}{Hemoglobin A1c (HbA1c) % (percentage).}
#'   \item{\code{Chol}}{Cholesterol (Chol) (mg/dL).}
#'   \item{\code{TG}}{Triglycerides (TG) (mg/dL).}
#'   \item{\code{HDL}}{High-density lipoprotein (HDL) (mg/dL).}
#'   \item{\code{LDL}}{Low-density lipoprotein (LDL) (mg/dL).}
#'   \item{\code{VLDL}}{Very-low-density lipoprotein (VLDL) (mg/dL).}
#'   \item{\code{BMI}}{Body mass index (BMI).}
#'   \item{\code{CLASS}}{Class (the patient's diabetes disease class may be Diabetic, Non-Diabetic, or Predict-Diabetic).}
#' }
#' @details
#' The data were collected from the Iraqi society, as they data were acquired from the laboratory of Medical City Hospital
#' and (the Specializes Center for Endocrinology and Diabetes-Al-Kindy Teaching Hospital).
#'
#' @source Rashid A (2020), “Diabetes Dataset”, Mendeley Data, V1, doi: 10.17632/wj9rwkp9c2.1.
#' @examples
#' data(diabetes)
#' knitr::kable(head(diabetes,5),caption="Five individuals in diabetes data")
#' @seealso \code{\link[gaawr2]{DiaHealth}}.

"diabetes"

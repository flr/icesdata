#' Master stock list for BIM / multi-fishery case studies.
#'
#' ICES FishStock codes (plus ICCAT \code{alb-n}) with fishery grouping and
#' short display names. Used as the default \code{sid} set for
#' \code{\link{getSAG}} when \code{sid = NULL}.
#'
#' @return Data frame with columns \code{fishery}, \code{sid}, and \code{common}.
#' @export
#' @examples
#' head(bimSids())
#' subset(bimSids(), fishery == "Pelagics")
bimSids <- function() {
  data.frame(
    fishery = c(
      rep("Pelagics", 4),
      "ICCAT",
      rep("Demersal", 8),
      rep("Nephrops", 9)
    ),
    sid = c(
      "boc.27.6-8",
      "hom.27.2a3a4a5b6a7a-ce-k8",
      "mac.27.nea",
      "whb.27.1-91214",
      "alb-n",
      "anf.27.3a46",
      "ank.27.78abd",
      "hke.27.3a46-8abd",
      "meg.27.7b-k8abd",
      "meg.27.8c9a",
      "mon.27.78abd",
      "whg.27.47d",
      "whg.27.7b-ce-k",
      "nep.fu.11",
      "nep.fu.12",
      "nep.fu.13",
      "nep.fu.14",
      "nep.fu.15",
      "nep.fu.16",
      "nep.fu.19",
      "nep.fu.2021",
      "nep.fu.22"
    ),
    common = c(
      "Boarfish",
      "Horse mackerel",
      "Mackerel",
      "Blue whiting",
      "Northern albacore",
      "Anglerfish (3a46)",
      "Black anglerfish",
      "Northern hake",
      "Megrim (7b-k8)",
      "Megrim (8c9a)",
      "White anglerfish",
      "Whiting (4+7d)",
      "Whiting (7b-ce-k)",
      "Nephrops FU 11",
      "Nephrops FU 12",
      "Nephrops FU 13",
      "Nephrops FU 14",
      "Nephrops FU 15",
      "Nephrops FU 16",
      "Nephrops FU 19",
      "Nephrops FU 20–21",
      "Nephrops FU 22"
    ),
    stringsAsFactors = FALSE
  )
}

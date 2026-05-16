# =============================================================================
# Cubic Spline Interpolation - Interactive Shiny App
# Numerical Methods Final Project
#
# Authors:      
#   - Francis Jay Abordo
#   - Clybel Djen Bonachita
#   - Jake Harvey Despabeladero
#   - Kent Anthony Dulangon
#   - Vince Quijano
# Date:         May 2026
#
# WHAT THIS APP DOES
#   Builds the natural cubic spline through a user-supplied set of (x, y)
#   data points, displays the resulting smooth curve, and estimates the
#   y-value at any chosen x. Three top-level tabs:
#     - Introduction : plain-language explanation of cubic splines, with
#                      side-by-side intro plots (straight lines vs. spline)
#     - Calculator   : enter points, see the interactive curve, watch it
#                      assemble one piece at a time, and walk through the
#                      seven-step derivation
#     - Solution     : MathJax-rendered piecewise polynomials, coefficient
#                      table, and the combined piecewise function
#
# METHOD
#   Natural cubic spline. For n+1 points (x_0, y_0) ... (x_n, y_n) we build
#   n cubic pieces
#       S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
#   such that adjacent pieces agree in value, first derivative, and second
#   derivative at every interior point, with the natural boundary conditions
#   c_0 = c_n = 0. The c_i are obtained by solving a tridiagonal linear
#   system; a_i, b_i, d_i then follow in closed form. The solver is written
#   out by hand in compute_spline() below (we deliberately do NOT call R's
#   built-in splinefun() - the manual implementation is the educational
#   point of the project and is exposed in the Calculator and Solution tabs).
#
# HOW TO RUN
#   1. Install dependencies (one-time):
#        install.packages(c("shiny", "plotly"))
#   2. From this directory:
#        Rscript -e 'shiny::runApp(".", launch.browser = TRUE)'
#      or open this file in RStudio and click "Run App".
#
# DEPENDENCIES
#   shiny, plotly. MathJax is loaded via shiny::withMathJax() - no extra
#   install needed.
#
# FILE ROADMAP  (approximate line numbers)
#     1-  50 : header + helpers + the natural-cubic-spline solver
#    50- 600 : UI definition (CSS theme, hero, three tabs)
#   600- end : server logic (reactives, plots, step-by-step derivation,
#              presentation-friendly Solution view)
# =============================================================================

library(shiny)
library(plotly)

# ---- Small UI helper: tooltip icon (uses native browser title attr) ----
help_tip <- function(text) {
  tags$span(class = "help-icon", title = text, HTML("&#9432;"))
}

# ============================================================================
# BEGIN experiment: SVG illustrations for "Where is it used?" bullets.
# If this looks out of place, delete this whole block AND its matching CSS
# block (search for "BEGIN experiment" in the stylesheet) AND change the
# tags$ul class back from "uses-list uses-icons" to just "uses-list" and
# remove the five use_case_icon(...) calls in the corresponding tags$li.
# ============================================================================
use_case_icon <- function(kind) {
  svg <- switch(kind,
    "animation" = paste0(
      "<svg viewBox='0 0 36 36' width='28' height='28' aria-hidden='true'>",
      "<path d='M 4 28 Q 14 4, 22 18 T 32 8' stroke='#3FB6A1' ",
        "stroke-width='2.5' fill='none' stroke-linecap='round'/>",
      "<circle cx='32' cy='8' r='3' fill='#F26D80'/>",
      "</svg>"),
    "typography" = paste0(
      "<svg viewBox='0 0 36 36' width='28' height='28' aria-hidden='true'>",
      "<text x='18' y='29' text-anchor='middle' font-family='Georgia, serif' ",
        "font-size='30' font-style='italic' font-weight='600' fill='#3FB6A1'>",
        "&amp;</text>",
      "</svg>"),
    "cad" = paste0(
      "<svg viewBox='0 0 36 36' width='28' height='28' aria-hidden='true'>",
      "<path d='M 3 24 Q 6 15, 12 14 L 22 14 Q 28 14, 33 24 Z' ",
        "fill='#3FB6A1'/>",
      "<circle cx='10' cy='27' r='3' fill='#2D3748'/>",
      "<circle cx='26' cy='27' r='3' fill='#2D3748'/>",
      "</svg>"),
    "data" = paste0(
      "<svg viewBox='0 0 36 36' width='28' height='28' aria-hidden='true'>",
      "<path d='M 4 26 Q 10 8, 18 20 T 32 10' stroke='#3FB6A1' ",
        "stroke-width='2.5' fill='none' stroke-linecap='round'/>",
      "<circle cx='4' cy='26' r='2.5' fill='#F26D80'/>",
      "<circle cx='18' cy='20' r='2.5' fill='#F26D80'/>",
      "<circle cx='32' cy='10' r='2.5' fill='#F26D80'/>",
      "</svg>"),
    "engineering" = paste0(
      "<svg viewBox='0 0 36 36' width='28' height='28' aria-hidden='true'>",
      "<path d='M 4 20 Q 8 6, 12 26 T 20 14 T 28 24 T 32 10' ",
        "stroke='#A0AEC0' stroke-width='1.5' fill='none' stroke-dasharray='3 2'/>",
      "<path d='M 4 20 Q 18 4, 32 10' stroke='#3FB6A1' ",
        "stroke-width='2.5' fill='none' stroke-linecap='round'/>",
      "</svg>"))
  tags$span(class = "use-icon", HTML(svg))
}
# ============================================================================
# END experiment: SVG illustrations.
# ============================================================================

# ---- Natural Cubic Spline (manual implementation) ----
# Solves the tridiagonal system for the second-derivative coefficients
# under natural boundary conditions (c_0 = c_n = 0), then derives
# a_i, b_i, c_i, d_i so each interval has
#   S_i(x) = a_i + b_i (x - x_i) + c_i (x - x_i)^2 + d_i (x - x_i)^3.
#
# Numerical correctness check
# ---------------------------
# This manual solver has been verified against R's built-in
# splinefun(method = "natural") on the test case
#     x <- c(1, 2, 3, 4, 5);  y <- c(1, 4, 2, 5, 3)
#
#   x_query   manual              splinefun           |diff|
#   1.0       1.0000000000        1.0000000000        0
#   1.5       3.1696428571        3.1696428571        0
#   2.0       4.0000000000        4.0000000000        0
#   2.5       2.8660714286        2.8660714286        4.4e-16
#   3.0       2.0000000000        2.0000000000        4.4e-16
#   3.5       3.3660714286        3.3660714286        4.4e-16
#   4.0       5.0000000000        5.0000000000        0
#   4.5       4.6696428571        4.6696428571        8.9e-16
#   5.0       3.0000000000        3.0000000000        4.4e-16
#
# Maximum |manual - splinefun| over [1, 5] sampled at 1001 points: 1.78e-15
# (i.e. the two implementations agree to floating-point precision).
compute_spline <- function(x, y) {
  n <- length(x) - 1
  h <- diff(x)

  if (n >= 2) {
    A   <- matrix(0, n - 1, n - 1)
    rhs <- numeric(n - 1)
    for (i in 1:(n - 1)) {
      if (i > 1)     A[i, i - 1] <- h[i]
      A[i, i]        <- 2 * (h[i] + h[i + 1])
      if (i < n - 1) A[i, i + 1] <- h[i + 1]
      rhs[i] <- 3 * ((y[i + 2] - y[i + 1]) / h[i + 1] -
                     (y[i + 1] - y[i])     / h[i])
    }
    c_inner <- solve(A, rhs)
  } else {
    A <- NULL; rhs <- NULL; c_inner <- numeric(0)
  }
  c_full <- c(0, c_inner, 0)

  a <- y[1:n]
  b <- numeric(n)
  d <- numeric(n)
  for (i in 1:n) {
    b[i] <- (y[i + 1] - y[i]) / h[i] -
            h[i] * (2 * c_full[i] + c_full[i + 1]) / 3
    d[i] <- (c_full[i + 1] - c_full[i]) / (3 * h[i])
  }

  list(x = x, y = y, h = h,
       a = a, b = b, c = c_full[1:n], d = d,
       c_full = c_full, A = A, rhs = rhs)
}

eval_spline <- function(sp, xq) {
  n <- length(sp$h)
  i <- max(1, min(n, findInterval(xq, sp$x, all.inside = TRUE)))
  dx <- xq - sp$x[i]
  sp$a[i] + sp$b[i] * dx + sp$c[i] * dx^2 + sp$d[i] * dx^3
}

# ============================================================================
# UI section builders
# Each create_*() function returns one chunk of the UI tree. The
# ui <- fluidPage(...) composition at the bottom of this section just stitches
# them together, so the overall page layout stays readable at a glance.
# ============================================================================

# ---- Optional hero image at the top of the Introduction tab ----
# Leave as NULL to skip; paste a base64 data URI here to enable. To convert
# an image to base64 from R:
#   base64enc::dataURI(file = "hero.jpg", mime = "image/jpeg")
# Recommended dimensions: ~1600x500 (or any wide landscape crop). Keep file
# size under ~100 KB after encoding to keep app.R reasonable.
hero_image_uri <- NULL  # e.g. "data:image/jpeg;base64,/9j/4AAQSk..."

create_styles <- function() {
  tags$style(HTML("
    /* ---- base palette ----
       coral  #F26D80   teal   #3FB6A1
       cream  #FAFAF7   ink    #2D3748
       muted  #718096   line   #E2E6EC                            */

    body, .container-fluid {
      background: #FAFAF7;
      color: #2D3748;
      font-family: 'Roboto', system-ui, -apple-system, sans-serif;
      font-weight: 400;
    }
    h1, h2, h3, h4, h5 { font-weight: 500; color: #2D3748; }
    h1 { font-size: 30px; }
    h2 { font-size: 24px; }
    h3 { font-size: 20px; }
    h4 { font-size: 17px; }
    p  { line-height: 1.55; }
    a  { color: #F26D80; }

    /* Inline emphasis for key terms in prose -- soft coral marker-pen wash
       under the text. Use with span(class = 'highlighted-text', 'term'). */
    .highlighted-text {
      background: linear-gradient(180deg,
                                  transparent 62%,
                                  rgba(242, 109, 128, 0.22) 62%);
      font-weight: 500;
      padding: 0 2px;
    }

    /* Title panel */
    .title-hero { padding: 28px 20px 6px; }
    .title-hero h1 { margin: 0; font-weight: 600; }
    .title-hero .tagline { color: #718096; margin-top: 6px; font-size: 16px; }

    /* Optional hero image at the top of the Introduction tab */
    .intro-hero-image {
      margin-bottom: 18px;
      border-radius: 12px;
      overflow: hidden;
      box-shadow: 0 2px 10px rgba(0, 0, 0, 0.06);
    }
    .intro-hero-image img {
      width: 100%;
      height: auto;
      display: block;
    }

    /* Tabs */
    .nav-tabs { border-bottom: 2px solid #E2E6EC; margin-bottom: 18px; }
    .nav-tabs > li > a {
      color: #718096; border: none; padding: 12px 22px;
      font-weight: 500; font-size: 15px;
      background: transparent !important; border-radius: 0;
    }
    .nav-tabs > li > a:hover { color: #2D3748; background: transparent; }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:hover,
    .nav-tabs > li.active > a:focus {
      color: #F26D80; border: none;
      box-shadow: inset 0 -3px 0 #F26D80;
      background: transparent !important;
    }

    /* Inner sub-tabs (Calculator > Result & Plot / Build / How it works)
       styled as pill buttons in a soft tray so they're obviously clickable. */
    .calc-subtabs-hint {
      color: #718096; font-size: 13px;
      letter-spacing: 0.04em; text-transform: uppercase;
      margin: 4px 4px 6px; font-weight: 600;
    }
    .calc-subtabs .nav-tabs {
      display: flex; gap: 6px;
      border-bottom: none;
      background: #F1F4F9;
      padding: 6px;
      border-radius: 12px;
      margin-bottom: 18px;
    }
    .calc-subtabs .nav-tabs > li { margin-bottom: 0; }
    .calc-subtabs .nav-tabs > li > a {
      border-radius: 8px;
      padding: 10px 18px;
      background: transparent !important;
      box-shadow: none !important;
      color: #4A5568;
    }
    .calc-subtabs .nav-tabs > li > a:hover {
      background: #FFFFFF !important;
      color: #2D3748;
    }
    .calc-subtabs .nav-tabs > li.active > a,
    .calc-subtabs .nav-tabs > li.active > a:hover,
    .calc-subtabs .nav-tabs > li.active > a:focus {
      background: #FFFFFF !important;
      color: #F26D80;
      box-shadow: 0 2px 6px rgba(0, 0, 0, 0.08) !important;
    }

    /* Cards */
    .card {
      background: #FFFFFF;
      border-radius: 12px;
      box-shadow: 0 2px 10px rgba(0, 0, 0, 0.06);
      padding: 22px 26px;
      margin-bottom: 18px;
    }
    .card h3:first-child, .card h4:first-child { margin-top: 0; }
    .card-hero { background: linear-gradient(135deg, #FFF6F7 0%, #F0FBF8 100%);
                 border: 1px solid #FCE0E5; }
    .card-hero h2 { margin: 0 0 4px; color: #2D3748; }
    .card-hero .tag { color: #718096; font-size: 16px; }

    /* Sidebar override -- treat it as a card too */
    .well {
      background: #FFFFFF !important;
      border: none !important;
      border-radius: 12px !important;
      box-shadow: 0 2px 10px rgba(0, 0, 0, 0.06) !important;
      padding: 22px !important;
    }

    /* Buttons */
    .btn-default {
      background: #FFFFFF; color: #2D3748;
      border: 1px solid #E2E6EC; border-radius: 6px;
      padding: 6px 14px; font-weight: 500;
    }
    .btn-default:hover { background: #F7F8FB; border-color: #CBD2DA; }
    .btn-primary-coral {
      background: #F26D80; color: #FFFFFF !important;
      border: none; border-radius: 8px;
      padding: 11px 24px; font-weight: 500; font-size: 15px;
      box-shadow: 0 2px 6px rgba(242, 109, 128, 0.35);
      transition: transform 0.08s ease, box-shadow 0.12s ease;
    }
    .btn-primary-coral:hover, .btn-primary-coral:focus {
      background: #EC5A6F; color: #FFFFFF;
      transform: translateY(-1px);
      box-shadow: 0 4px 10px rgba(242, 109, 128, 0.45);
    }

    /* Soft teal action button (e.g. Load example) */
    .btn-soft-teal {
      background: #E8F6F1; color: #1F6B5B;
      border: 1px solid #BDE2D7; border-radius: 8px;
      padding: 8px 14px; font-weight: 500; font-size: 14px;
    }
    .btn-soft-teal:hover, .btn-soft-teal:focus {
      background: #D7EFE6; color: #1F6B5B;
      border-color: #A0D9C6;
    }

    /* Import / export widgets */
    .csv-hint {
      margin-top: -8px; margin-bottom: 12px;
      color: #718096; font-size: 12px; line-height: 1.4;
    }
    .csv-hint code {
      background: #F7F8FB; color: #2D3748;
      padding: 1px 5px; border-radius: 3px; font-size: 11px;
    }
    .btn-block-download {
      display: block; width: 100%;
      text-align: center;
    }

    /* Tooltip icon */
    .help-icon {
      display: inline-block; width: 18px; height: 18px;
      line-height: 18px; text-align: center;
      background: #E5F2EE; color: #3FB6A1;
      border-radius: 50%; font-size: 12px;
      cursor: help; margin-left: 4px;
      vertical-align: middle;
    }
    .help-icon:hover { background: #3FB6A1; color: #FFFFFF; }

    /* Data points table */
    .points-table .form-group { margin-bottom: 4px; }
    .points-table input { padding: 4px 6px; height: 30px;
                           border-radius: 4px; border: 1px solid #E2E6EC; }
    .points-table input:focus { border-color: #3FB6A1;
                                 box-shadow: 0 0 0 2px rgba(63, 182, 161, 0.18); }
    .points-table .row-num { padding-top: 6px; color: #A0AEC0;
                              font-size: 13px; text-align: center; }
    .points-table .header { font-weight: 600; padding: 4px 0;
                            color: #718096; font-size: 13px;
                            text-transform: uppercase; letter-spacing: 0.04em;
                            border-bottom: 1px solid #E2E6EC;
                            margin-bottom: 6px; text-align: center; }

    /* Step boxes */
    .step-box {
      border: 1px solid #E2E6EC; border-radius: 10px;
      padding: 14px 20px; margin-bottom: 14px;
      background: #FFFFFF;
      box-shadow: 0 1px 4px rgba(0, 0, 0, 0.04);
    }
    .step-box h4 { margin-top: 4px; color: #2D3748; }
    .step-box .step-friendly {
      color: #718096; font-size: 14px;
      margin: -2px 0 8px; font-style: italic;
    }
    .step-box pre {
      background: #F7F8FB; border: 1px solid #E2E6EC;
      border-radius: 6px; font-size: 13px; padding: 10px 12px;
    }
    .step-box.solution-step {
      background: #FFFBEC; border-color: #F2D27A;
    }
    .step-box.solution-step h4 { color: #7A5A00; }
    /* Collapsible step boxes (each <details> opens by default; users can
       click the heading to collapse a step they have already understood). */
    details.step-box > summary {
      list-style: none;
      cursor: pointer;
      padding-right: 28px;
      position: relative;
    }
    details.step-box > summary::-webkit-details-marker { display: none; }
    details.step-box > summary::after {
      content: '▾';  /* down-pointing arrow */
      position: absolute;
      right: 4px; top: 10px;
      color: #A0AEC0;
      transition: transform 0.2s ease;
      font-size: 14px;
    }
    details.step-box[open] > summary::after { transform: rotate(180deg); }
    details.step-box > summary h4 {
      display: inline-block;
      margin: 0;
      vertical-align: middle;
    }
    details.step-box > summary:hover h4 { color: #F26D80; }
    details.step-box.solution-step > summary:hover h4 { color: #7A5A00; }

    /* Solution sections + coefficient table */
    .solution-section { margin-bottom: 24px; }
    .solution-section h4 {
      color: #2D3748; padding-bottom: 6px;
      box-shadow: inset 0 -2px 0 #F26D80;
      display: inline-block; padding-right: 8px;
    }
    .coef-table { width: 100%; border-collapse: collapse;
                  margin: 8px 0 12px; }
    .coef-table th, .coef-table td {
      border: 1px solid #E2E6EC;
      padding: 8px 12px; text-align: right;
      font-family: 'JetBrains Mono', 'Courier New', monospace;
    }
    .coef-table th {
      background: #FCEEF1; font-family: 'Roboto', sans-serif;
      color: #7A2A3A; text-align: center; font-weight: 500;
    }
    .coef-table tr:nth-child(even) td { background: #FAFAF7; }

    /* Result / extrapolation callouts */
    .solution-callout {
      background: #E8F6F1; border: 1px solid #3FB6A1;
      border-left: 5px solid #3FB6A1; border-radius: 8px;
      padding: 16px 22px; margin: 16px 0; font-size: 15px;
    }
    .solution-callout strong { color: #1F6B5B; }
    .result-callout {
      background: #FFF1F3; border: 1px solid #F26D80;
      border-left: 5px solid #F26D80; border-radius: 8px;
      padding: 16px 22px; margin: 8px 0; font-size: 18px;
      font-weight: 500; color: #7A2A3A;
    }

    /* Calculation breakdown table -- ledger style, not code-y */
    .calc-table {
      width: 100%; border-collapse: collapse;
      margin: 6px 0 2px;
      font-family: 'Roboto', sans-serif;
      font-size: 14.5px;
    }
    .calc-table th, .calc-table td {
      padding: 11px 14px;
      border-bottom: 1px solid #EEF0F4;
      vertical-align: middle;
    }
    .calc-table thead th {
      background: transparent;
      color: #A0AEC0;
      font-weight: 600;
      font-size: 11px;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      text-align: left;
      border-bottom: 2px solid #E2E6EC;
      padding-bottom: 8px;
    }
    .calc-table thead th.value-col { text-align: right; }
    .calc-table td.step-num {
      width: 28px; color: #CBD5E0; font-weight: 600;
      text-align: center;
    }
    .calc-table td.step-label {
      color: #2D3748; font-weight: 500; width: 38%;
    }
    .calc-table td.step-expr {
      color: #718096;
    }
    .calc-table td.step-expr sub,
    .calc-table td.step-expr sup { font-size: 75%; }
    .calc-table td.step-value {
      text-align: right;
      font-variant-numeric: tabular-nums;
      color: #2D3748; font-weight: 500;
      white-space: nowrap;
    }
    .calc-table td.step-value.pos { color: #2D6A4F; }
    .calc-table td.step-value.neg { color: #B0455C; }
    .calc-table tr.total td {
      border-top: 2px solid #F26D80;
      border-bottom: none;
      padding-top: 14px;
    }
    .calc-table tr.total td.step-label { color: #2D3748; font-weight: 600; }
    .calc-table tr.total td.step-value {
      color: #F26D80; font-weight: 700; font-size: 18px;
    }
    .calc-empty {
      color: #A0AEC0; font-style: italic;
      padding: 14px 4px 4px;
    }

    /* Build-piece view */
    .calc-table td.interval-cell {
      color: #2D3748; font-weight: 500;
      font-variant-numeric: tabular-nums;
      white-space: nowrap;
    }
    .calc-table tr.latest-piece td {
      background: #FFF1F3;
    }
    .calc-table tr.latest-piece td.step-num {
      color: #F26D80; font-weight: 700;
    }
    .calc-table tr.latest-piece td.interval-cell {
      color: #7A2A3A; font-weight: 600;
    }

    .build-current {
      margin-top: 10px;
      padding: 14px 18px;
      background: #FFF6F7;
      border: 1px solid #FCE0E5;
      border-left: 4px solid #F26D80;
      border-radius: 8px;
      line-height: 1.6;
    }
    .build-current .formula-label {
      display: block; margin-top: 8px;
      color: #A0AEC0; font-size: 11px;
      letter-spacing: 0.08em; text-transform: uppercase;
    }
    .build-current .formula {
      font-size: 16px; color: #2D3748;
      margin-top: 2px;
    }
    .build-empty {
      color: #A0AEC0; font-style: italic;
      padding: 14px 0;
    }
    .point-dot {
      display: inline-block; width: 9px; height: 9px;
      background: #F26D80; border-radius: 50%;
      margin-right: 4px; vertical-align: middle;
    }

    /* Caption under intro mini-plots */
    .plot-caption {
      text-align: center; color: #718096;
      font-size: 13px; margin: 4px 0 12px;
    }

    /* Numbered how-to list */
    .howto-list { counter-reset: how; padding-left: 0; list-style: none; }
    .howto-list li {
      counter-increment: how;
      position: relative; padding: 4px 0 6px 38px;
      margin-bottom: 4px;
    }
    .howto-list li::before {
      content: counter(how); position: absolute; left: 0; top: 2px;
      width: 26px; height: 26px; line-height: 26px;
      background: #3FB6A1; color: #FFFFFF;
      border-radius: 50%; text-align: center; font-weight: 600;
      font-size: 13px;
    }

    /* Use-case bullets with dot */
    .uses-list { padding-left: 0; list-style: none; }
    .uses-list li {
      padding: 6px 0 6px 20px; position: relative;
    }
    .uses-list li::before {
      content: ''; position: absolute; left: 0; top: 14px;
      width: 8px; height: 8px; border-radius: 50%;
      background: #F26D80;
    }
    /* BEGIN experiment: SVG illustrations for use-case bullets.
       Remove this block (down to the END experiment marker) to fully revert. */
    .uses-list.uses-icons li {
      padding: 8px 0 8px 44px; min-height: 36px;
    }
    .uses-list.uses-icons li::before { display: none; }
    .uses-list.uses-icons .use-icon {
      position: absolute; left: 0; top: 4px;
      width: 32px; height: 32px;
      display: inline-flex; align-items: center; justify-content: center;
    }
    /* END experiment */
  "))
}

create_hero <- function() {
  div(class = "title-hero",
    h1("Cubic Spline Interpolation"),
    div(class = "tagline",
      "Draw a smooth curve through any set of points — ",
      "an interactive guide to a classic numerical method.")
  )
}

create_intro_tab <- function() {
  tabPanel("Introduction",
      fluidRow(
        column(width = 10, offset = 1,
          # ---- Optional hero image (renders only if hero_image_uri is set) ----
          if (!is.null(hero_image_uri)) {
            div(class = "intro-hero-image",
              img(src = hero_image_uri,
                  alt = "Cubic spline interpolation in the real world"))
          },
          div(class = "card card-hero",
            h2("What is a cubic spline?"),
            div(class = "tag",
              "Imagine connecting dots on a graph — but instead of straight ",
              "lines that kink at every point, you want one perfectly smooth ",
              "curve. That's what a ",
              span(class = "highlighted-text", "cubic spline"),
              " does."
            )
          ),

          div(class = "card",
            h3("See it in action"),
            p(style = "color:#718096;",
              "Same five points, two different ways of joining them. ",
              "Look at the difference:"),
            fluidRow(
              column(6,
                plotlyOutput("intro_plot_lines", height = "240px"),
                div(class = "plot-caption",
                  HTML("<strong>Straight lines</strong> &mdash; sharp corners at every point.")
                )
              ),
              column(6,
                plotlyOutput("intro_plot_spline", height = "240px"),
                div(class = "plot-caption",
                  HTML("<strong>Cubic spline</strong> &mdash; smooth curves the whole way.")
                )
              )
            )
          ),

          div(class = "card",
            h3("How does it work?"),
            p(withMathJax(HTML(paste0(
              "Between every pair of your data points, the spline lays down ",
              "a small cubic polynomial — a curve of the form ",
              "\\(a + bx + cx^2 + dx^3\\). ",
              "Each piece is chosen so that where two pieces meet, they share ",
              "<span class='highlighted-text'>the same value, the same slope, ",
              "and the same curvature</span>. ",
              "That's why the joins are invisible.")))),
            p(style = "color:#718096; font-size: 14px;",
              "The flavour used here is the ",
              span(class = "highlighted-text", "natural"),
              " cubic spline: at the two endpoints the curvature is set to ",
              "zero, so the curve tapers off gently rather than flicking up."),
            h3("Why does it work?"),
            p("Many real-world relationships are not linear. Cubic splines ",
              "work because they enforce ",
              span(class = "highlighted-text", "smooth connections between cubic segments"),
              ", which naturally produces accurate approximations ",
              "of real functions and minimises errors that sharp or straight ",
              "lines introduce. That makes downstream calculations like ",
              "area, slope, and motion much more accurate and stable.")
          ),

          div(class = "card",
            h3("Where is it used?"),
            tags$ul(class = "uses-list uses-icons",  # remove " uses-icons" + the use_case_icon() calls below to revert
              tags$li(use_case_icon("animation"),
                tags$strong("Animation & games — "),
                "smooth motion paths for characters and cameras."),
              tags$li(use_case_icon("typography"),
                tags$strong("Typography — "),
                "the curves of letterforms in modern fonts."),
              tags$li(use_case_icon("cad"),
                tags$strong("Computer-aided design — "),
                "shaping car bodies, aircraft wings, product surfaces."),
              tags$li(use_case_icon("data"),
                tags$strong("Data analysis — "),
                "filling in missing measurements between known samples."),
              tags$li(use_case_icon("engineering"),
                tags$strong("Engineering — "),
                "approximating complex functions with something easy to compute.")
            )
          ),

          div(class = "card",
            h3("Using this app"),
            tags$ol(class = "howto-list",
              tags$li(tags$strong("Type in your points."),
                " On the Calculator tab, fill the X / Y table with the ",
                "coordinates you want the curve to pass through."),
              tags$li(tags$strong("Pick an x to estimate."),
                " Tell the app which x-value you want a smooth-curve ",
                "y-value for."),
              tags$li(tags$strong("Read the answer."),
                " The estimate, an interactive plot, the full derivation, ",
                "and the polynomial formulas are all there for you.")
            ),
            div(style = "text-align:center; margin-top: 20px;",
              actionButton("go_to_calc",
                HTML("Try it yourself &rarr;"),
                class = "btn-primary-coral")
            )
          )
        )
      )
    )
}

create_calculator_tab <- function() {
    tabPanel("Calculator",
      sidebarLayout(
        sidebarPanel(
          h4("Your data points ",
             help_tip(paste0(
               "Each row is one (x, y) pair. The spline will pass through ",
               "every point you list here. Add or remove rows with the ",
               "buttons below."))),
          fluidRow(
            column(6, actionButton("add_row",    "+ Add row",    width = "100%")),
            column(6, actionButton("remove_row", "- Remove row", width = "100%"))
          ),
          br(),
          actionButton("load_example", "Load example data",
                       width = "100%", class = "btn-soft-teal"),
          br(), br(),
          fileInput("import_csv", "Import from CSV",
                    accept = c(".csv"),
                    buttonLabel = "Choose CSV",
                    placeholder = "No file selected"),
          tags$div(class = "csv-hint",
                   "Two columns, named ", tags$code("x"), " and ", tags$code("y"),
                   ", or any two numeric columns."),
          downloadButton("export_csv", "Download current points as CSV",
                         class = "btn-default btn-block-download"),
          br(), br(),
          uiOutput("point_inputs"),
          br(),
          tags$label("Estimate y at this x ",
                     help_tip(paste0(
                       "The app reads the smooth curve at this x value and ",
                       "tells you what y is there. Pick any x between your ",
                       "smallest and largest data point."))),
          numericInput("x_query", label = NULL, value = 3.5),
          helpText("You need at least 3 points, and every x-value must be unique.")
        ),

        mainPanel(
          div(class = "calc-subtabs-hint", "Switch view →"),
          div(class = "calc-subtabs",
          tabsetPanel(
            tabPanel("Result & Plot",
              div(class = "card",
                h4("Estimated value"),
                uiOutput("result_callout")
              ),
              div(class = "card",
                h4("Curve through your points"),
                p(style = "color:#718096; font-size:14px; margin-top:-4px;",
                  HTML("<strong>Tip:</strong> hover the curve to read values, ",
                       "drag a box to zoom in, double-click to reset.")),
                plotlyOutput("splinePlot", height = "450px")
              )
            ),
            tabPanel("Build piece by piece",
              br(),
              div(class = "card",
                h4("Watch the spline assemble"),
                p(style = "color:#718096; font-size:14px; margin-top:-4px;",
                  HTML("A cubic spline isn't built by iterating until something ",
                       "converges &mdash; it's built piece by piece, one cubic ",
                       "polynomial per gap between your points. Move the slider ",
                       "below to add one piece at a time and watch the curve ",
                       "grow.")),
                sliderInput("piece_count",
                            label = "Pieces drawn:",
                            min = 1, max = 6, value = 1, step = 1,
                            width = "100%"),
                plotlyOutput("build_plot", height = "400px")
              ),
              div(class = "card",
                h4("Latest piece"),
                p(style = "color:#718096; font-size:14px; margin-top:-4px;",
                  "The cubic polynomial that fills the most recently added gap."),
                uiOutput("build_current")
              ),
              div(class = "card",
                h4("Pieces built so far"),
                p(style = "color:#718096; font-size:14px; margin-top:-4px;",
                  "Each row is one cubic piece of the spline, with its width ",
                  HTML("(<em>h</em>) and the four coefficients ",
                       "<em>a, b, c, d</em>.")),
                uiOutput("build_table")
              )
            ),
            tabPanel("How it works",
              br(),
              uiOutput("steps_ui")
            )
          )
          )
        )
      )
    )
}

create_solution_tab <- function() {
    tabPanel("Solution",
      fluidRow(
        column(width = 10, offset = 1,
          withMathJax(),
          br(),
          uiOutput("solution_ui")
        )
      )
    )
}

# ---- UI composition ----
ui <- fluidPage(
  withMathJax(),
  tags$head(
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;600;700&display=swap")
  ),
  create_styles(),
  create_hero(),
  tabsetPanel(id = "main_tabs",
    create_intro_tab(),
    create_calculator_tab(),
    create_solution_tab()
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # Sample data used by the "Load example data" button (NOT shown on first load).
  default_x <- c(1, 2, 3, 4, 5, 6, 7)
  default_y <- c(2, 3, 5, 4, 6, 8, 7)
  # The app opens with 7 blank rows -- no pre-filled numbers, so visitors are
  # not confused into thinking the example output is "their" answer.
  n_points  <- reactiveVal(7)

  observeEvent(input$add_row,    n_points(n_points() + 1))
  observeEvent(input$remove_row, {
    if (n_points() > 3) n_points(n_points() - 1)
  })

  # Helper: push a vector of x/y values into the visible numericInputs.
  # Used by both "Load example" and the CSV import handler so that the
  # values reliably overwrite anything already in the table.
  push_points <- function(x, y) {
    stopifnot(length(x) == length(y))
    n_points(length(x))
    for (i in seq_along(x)) {
      updateNumericInput(session, sprintf("x_%d", i), value = x[i])
      updateNumericInput(session, sprintf("y_%d", i), value = y[i])
    }
  }

  # One-click: populate the table with the sample data.
  observeEvent(input$load_example, push_points(default_x, default_y))

  # CSV import: read a two-column CSV (x, y) and populate the inputs.
  observeEvent(input$import_csv, {
    file <- input$import_csv
    req(file)
    df <- tryCatch(read.csv(file$datapath, header = TRUE,
                            stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) {
      showNotification("Unable to read CSV file.", type = "error")
      return()
    }

    # Accept files with columns named x/y (case-insensitive) or, failing
    # that, fall back to the first two columns.
    names_lower <- tolower(names(df))
    if (all(c("x", "y") %in% names_lower)) {
      xi <- which(names_lower == "x")[1]
      yi <- which(names_lower == "y")[1]
      pts <- data.frame(x = suppressWarnings(as.numeric(df[[xi]])),
                        y = suppressWarnings(as.numeric(df[[yi]])))
    } else if (ncol(df) >= 2) {
      pts <- data.frame(x = suppressWarnings(as.numeric(df[[1]])),
                        y = suppressWarnings(as.numeric(df[[2]])))
    } else {
      showNotification("CSV must have at least two numeric columns.",
                       type = "error")
      return()
    }

    if (any(is.na(pts$x)) || any(is.na(pts$y))) {
      showNotification("CSV contains non-numeric values.", type = "error")
      return()
    }
    if (nrow(pts) < 3) {
      showNotification("CSV needs at least 3 points.", type = "error")
      return()
    }
    if (any(duplicated(pts$x))) {
      showNotification("CSV contains duplicate x values.", type = "error")
      return()
    }

    pts <- pts[order(pts$x), , drop = FALSE]
    push_points(pts$x, pts$y)
    showNotification(sprintf("Imported %d points.", nrow(pts)),
                     type = "message")
  })

  # Render the X/Y input table. Each cell starts blank (NA); user typing or
  # one of the loaders (Load example / CSV import) writes values in.
  output$point_inputs <- renderUI({
    n <- n_points()
    isolate({
      get_val <- function(i, prefix) {
        v <- input[[paste0(prefix, "_", i)]]
        if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v))
      }
      cur_x <- vapply(seq_len(n), get_val, numeric(1), prefix = "x")
      cur_y <- vapply(seq_len(n), get_val, numeric(1), prefix = "y")
    })

    header <- fluidRow(class = "header",
      column(2, "#"),
      column(5, "X"),
      column(5, "Y")
    )
    rows <- lapply(seq_len(n), function(i) {
      fluidRow(
        column(2, div(class = "row-num", i)),
        column(5, numericInput(paste0("x_", i), NULL, value = cur_x[i])),
        column(5, numericInput(paste0("y_", i), NULL, value = cur_y[i]))
      )
    })
    div(class = "points-table", header, rows)
  })

  # Reactive: gather and validate user-entered points.
  data_points <- reactive({
    n <- n_points()
    read_col <- function(prefix) {
      vapply(seq_len(n), function(i) {
        v <- input[[paste0(prefix, "_", i)]]
        if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v))
      }, numeric(1))
    }
    x <- read_col("x")
    y <- read_col("y")

    validate(
      need(n >= 3, "Add at least 3 points to see the spline."),
      need(all(!is.na(x)), "Fill in every X value to see the spline."),
      need(all(!is.na(y)), "Fill in every Y value to see the spline."),
      need(!any(duplicated(x)), "X values must be unique.")
    )

    ord <- order(x)
    list(x = x[ord], y = y[ord])
  })

  # Reactive: build the natural cubic spline (manual solver)
  spline_data <- reactive({
    d <- data_points()
    compute_spline(d$x, d$y)
  })

  # Reactive: a vectorised evaluator wrapping the manual spline
  spline_fn <- reactive({
    sp <- spline_data()
    function(xq) sapply(xq, function(x) eval_spline(sp, x))
  })

  # Output: interpolated value as a coral callout
  output$result_callout <- renderUI({
    d <- data_points()
    f <- spline_fn()
    xq <- input$x_query

    validate(need(!is.na(xq), "Please enter a valid x to estimate."))

    yq <- f(xq)
    fmt <- function(v) formatC(round(v, 5), format = "g", digits = 5)

    if (xq < min(d$x) || xq > max(d$x)) {
      div(class = "result-callout",
          style = "background:#FFF1F3; color:#7A2A3A;",
          HTML(sprintf(
            "x = %s is outside your data range [%s, %s]. ",
            fmt(xq), fmt(min(d$x)), fmt(max(d$x)))),
          tags$div(style = "font-size:14px; font-weight:400; color:#A55067;",
                   "(That counts as extrapolation — the natural spline isn't reliable here.)"))
    } else {
      div(class = "result-callout",
          HTML(sprintf("f(%s) = <span style='color:#F26D80;'>%s</span>",
                       fmt(xq), fmt(yq))))
    }
  })

  output$export_csv <- downloadHandler(
    filename = function() sprintf("spline_points_%s.csv",
                                 format(Sys.time(), "%Y%m%d-%H%M%S")),
    content = function(file) {
      d <- data_points()
      df <- data.frame(x = d$x, y = d$y)
      write.csv(df, file, row.names = FALSE)
    }
  )

  # Output: interactive spline plot (plotly)
  output$splinePlot <- renderPlotly({
    d <- data_points()
    f <- spline_fn()
    xq <- input$x_query

    x_seq <- seq(min(d$x), max(d$x), length.out = 500)
    y_seq <- f(x_seq)

    p <- plot_ly() %>%
      add_trace(x = x_seq, y = y_seq, type = "scatter", mode = "lines",
                line = list(color = "#3FB6A1", width = 3),
                name = "Cubic spline",
                hovertemplate = "x = %{x:.4g}<br>y = %{y:.4g}<extra></extra>") %>%
      add_trace(x = d$x, y = d$y, type = "scatter", mode = "markers",
                marker = list(color = "#F26D80", size = 11,
                              line = list(color = "#FFFFFF", width = 2)),
                name = "Your points",
                hovertemplate = "Point<br>x = %{x}<br>y = %{y}<extra></extra>")

    if (!is.na(xq)) {
      yq <- tryCatch(f(xq), error = function(e) NA_real_)
      if (!is.na(yq) && xq >= min(d$x) && xq <= max(d$x)) {
        p <- p %>% add_trace(
          x = xq, y = yq, type = "scatter", mode = "markers",
          marker = list(color = "#2D6A4F", size = 14, symbol = "x",
                        line = list(color = "#2D6A4F", width = 3)),
          name = "Estimate",
          hovertemplate = "Estimate<br>x = %{x}<br>y = %{y:.4g}<extra></extra>")
      }
    }

    p %>% layout(
      xaxis = list(title = "x", gridcolor = "#EEF0F4",
                    zeroline = FALSE, showline = FALSE),
      yaxis = list(title = "y", gridcolor = "#EEF0F4",
                    zeroline = FALSE, showline = FALSE),
      paper_bgcolor = "#FFFFFF",
      plot_bgcolor  = "#FFFFFF",
      margin = list(l = 50, r = 20, t = 20, b = 50),
      legend = list(orientation = "h", y = -0.18, x = 0.5,
                    xanchor = "center", bgcolor = "rgba(0,0,0,0)"),
      font = list(family = "Roboto, sans-serif", color = "#2D3748")
    ) %>%
      config(displayModeBar = FALSE)
  })

  # ---- Build piece by piece ----

  # Keep the slider's max in sync with the number of intervals
  observe({
    d <- data_points()
    n <- length(d$x) - 1
    current <- isolate(input$piece_count)
    new_val <- if (is.null(current)) 1 else max(1, min(current, n))
    updateSliderInput(session, "piece_count",
                      min = 1, max = n, value = new_val)
  })

  output$build_plot <- renderPlotly({
    d  <- data_points()
    sp <- spline_data()
    n  <- length(sp$h)
    k  <- input$piece_count
    if (is.null(k)) k <- 1
    k  <- max(1, min(k, n))

    # Pre-compute full y range so plot doesn't rescale as pieces are added
    x_full <- seq(min(d$x), max(d$x), length.out = 300)
    y_full <- sapply(x_full, function(x) eval_spline(sp, x))
    y_pad  <- 0.08 * diff(range(c(d$y, y_full)))
    y_lo   <- min(d$y, y_full) - y_pad
    y_hi   <- max(d$y, y_full) + y_pad
    x_pad  <- 0.04 * diff(range(d$x))

    p <- plot_ly()

    # Built pieces: older ones teal, newest one coral & thicker
    for (i in seq_len(k)) {
      x_seg <- seq(sp$x[i], sp$x[i + 1], length.out = 80)
      y_seg <- sapply(x_seg, function(x) eval_spline(sp, x))
      is_latest <- (i == k)
      p <- add_trace(p,
        x = x_seg, y = y_seg,
        type = "scatter", mode = "lines",
        line = list(
          color = if (is_latest) "#F26D80" else "#3FB6A1",
          width = if (is_latest) 4.5 else 3
        ),
        showlegend = FALSE,
        hovertemplate = sprintf(
          "Piece S<sub>%d</sub><br>x = %%{x:.4g}<br>y = %%{y:.4g}<extra></extra>",
          i - 1))
    }

    # All data points
    p <- add_trace(p,
      x = d$x, y = d$y,
      type = "scatter", mode = "markers",
      marker = list(color = "#F26D80", size = 11,
                    line = list(color = "#FFFFFF", width = 2)),
      showlegend = FALSE,
      hovertemplate = "Point<br>x = %{x}<br>y = %{y}<extra></extra>")

    # Halo on the two endpoints of the latest piece
    p <- add_trace(p,
      x = c(sp$x[k], sp$x[k + 1]),
      y = c(sp$y[k], sp$y[k + 1]),
      type = "scatter", mode = "markers",
      marker = list(color = "rgba(242,109,128,0)", size = 22,
                    line = list(color = "#F26D80", width = 2)),
      showlegend = FALSE,
      hoverinfo = "skip")

    p %>% layout(
      xaxis = list(title = "x", gridcolor = "#EEF0F4",
                    zeroline = FALSE, showline = FALSE,
                    range = c(min(d$x) - x_pad, max(d$x) + x_pad)),
      yaxis = list(title = "y", gridcolor = "#EEF0F4",
                    zeroline = FALSE, showline = FALSE,
                    range = c(y_lo, y_hi)),
      paper_bgcolor = "#FFFFFF", plot_bgcolor = "#FFFFFF",
      margin = list(l = 50, r = 20, t = 20, b = 50),
      font = list(family = "Roboto, sans-serif", color = "#2D3748")
    ) %>%
      config(displayModeBar = FALSE)
  })

  output$build_current <- renderUI({
    sp <- spline_data()
    n  <- length(sp$h)
    k  <- input$piece_count
    if (is.null(k)) k <- 1
    k  <- max(1, min(k, n))
    i  <- k - 1

    fmt <- function(v) formatC(round(v, 5), format = "g", digits = 5)

    signed_term <- function(coef, expr) {
      if (abs(round(coef, 8)) < 1e-10) return("")
      if (coef >= 0) paste0(" + ", fmt(coef), expr)
      else            paste0(" &minus; ", fmt(-coef), expr)
    }

    xi    <- sp$x[k]
    xterm <- sprintf("(x &minus; %s)", fmt(xi))
    poly  <- paste0(
      fmt(sp$a[k]),
      signed_term(sp$b[k], xterm),
      signed_term(sp$c[k], paste0(xterm, "<sup>2</sup>")),
      signed_term(sp$d[k], paste0(xterm, "<sup>3</sup>"))
    )

    HTML(sprintf(paste0(
      "<div class='build-current'>",
      "<div>Piece <strong>%d of %d</strong> &middot; ",
      "called <strong>S<sub>%d</sub></strong> &middot; ",
      "connects <span class='point-dot'></span>(%s, %s) ",
      "to <span class='point-dot'></span>(%s, %s).</div>",
      "<span class='formula-label'>The cubic for this piece</span>",
      "<div class='formula'>S<sub>%d</sub>(x) = %s</div>",
      "</div>"),
      k, n,
      i,
      fmt(sp$x[k]),     fmt(sp$y[k]),
      fmt(sp$x[k + 1]), fmt(sp$y[k + 1]),
      i, poly))
  })

  output$build_table <- renderUI({
    sp <- spline_data()
    n  <- length(sp$h)
    k  <- input$piece_count
    if (is.null(k)) k <- 1
    k  <- max(1, min(k, n))

    fmt <- function(v) formatC(round(v, 5), format = "g", digits = 5)

    rows <- paste(sapply(seq_len(k), function(i) {
      cls <- if (i == k) " class='latest-piece'" else ""
      sprintf(paste0(
        "<tr%s>",
        "<td class='step-num'>%d</td>",
        "<td class='interval-cell'>[%s, %s]</td>",
        "<td class='step-value'>%s</td>",
        "<td class='step-value'>%s</td>",
        "<td class='step-value'>%s</td>",
        "<td class='step-value'>%s</td>",
        "<td class='step-value'>%s</td>",
        "</tr>"),
        cls, i - 1,
        fmt(sp$x[i]), fmt(sp$x[i + 1]),
        fmt(sp$h[i]),
        fmt(sp$a[i]), fmt(sp$b[i]), fmt(sp$c[i]), fmt(sp$d[i]))
    }), collapse = "")

    HTML(paste0(
      "<table class='calc-table'>",
      "<thead><tr>",
      "<th></th><th>Interval [x<sub>i</sub>, x<sub>i+1</sub>]</th>",
      "<th class='value-col'>Width <em>h</em></th>",
      "<th class='value-col'><em>a</em></th>",
      "<th class='value-col'><em>b</em></th>",
      "<th class='value-col'><em>c</em></th>",
      "<th class='value-col'><em>d</em></th>",
      "</tr></thead>",
      "<tbody>", rows, "</tbody></table>"
    ))
  })

  # ---- Intro tab: side-by-side comparison plots ----
  intro_x <- c(1, 2, 3, 4, 5)
  intro_y <- c(1, 4, 2, 5, 3)

  apply_intro_layout <- function(p) {
    layout(p,
      xaxis = list(title = "x", gridcolor = "#EEF0F4",
                    zeroline = FALSE, showline = FALSE, fixedrange = TRUE),
      yaxis = list(title = "y", gridcolor = "#EEF0F4",
                    zeroline = FALSE, showline = FALSE, fixedrange = TRUE),
      paper_bgcolor = "#FFFFFF",
      plot_bgcolor  = "#FFFFFF",
      margin = list(l = 40, r = 10, t = 10, b = 40),
      showlegend = FALSE,
      font = list(family = "Roboto, sans-serif", color = "#2D3748")
    )
  }

  output$intro_plot_lines <- renderPlotly({
    p <- plot_ly(x = intro_x, y = intro_y, type = "scatter",
                 mode = "lines+markers",
                 line   = list(color = "#F26D80", width = 3),
                 marker = list(color = "#F26D80", size = 11,
                               line = list(color = "#FFFFFF", width = 2)),
                 hoverinfo = "skip")
    config(apply_intro_layout(p), displayModeBar = FALSE)
  })

  output$intro_plot_spline <- renderPlotly({
    sp <- compute_spline(intro_x, intro_y)
    x_seq <- seq(min(intro_x), max(intro_x), length.out = 300)
    y_seq <- sapply(x_seq, function(x) eval_spline(sp, x))
    p <- plot_ly()
    p <- add_trace(p, x = x_seq, y = y_seq, type = "scatter", mode = "lines",
                   line = list(color = "#3FB6A1", width = 3),
                   hoverinfo = "skip")
    p <- add_trace(p, x = intro_x, y = intro_y, type = "scatter",
                   mode = "markers",
                   marker = list(color = "#F26D80", size = 11,
                                 line = list(color = "#FFFFFF", width = 2)),
                   hoverinfo = "skip")
    config(apply_intro_layout(p), displayModeBar = FALSE)
  })

  # Introduction / Learn tab CTA -> Calculator
  observeEvent(input$go_to_calc, {
    updateTabsetPanel(session, "main_tabs", selected = "Calculator")
  })

  # Output: step-by-step derivation of the spline
  output$steps_ui <- renderUI({
    sp <- spline_data()
    xq <- input$x_query
    n  <- length(sp$h)
    fmt <- function(v, k = 5) formatC(round(v, k), format = "fg", flag = "#")

    # Helper: build a small HTML <table> styled like a textbook coefficient
    # table. `headers` is a character vector of <th> contents (HTML allowed);
    # `body` is a list of character vectors, one per row, of <td> contents.
    mini_table <- function(headers, body, max_width = NULL) {
      header_html <- paste0("<tr>",
        paste(sprintf("<th>%s</th>", headers), collapse = ""),
        "</tr>")
      body_html <- paste(
        vapply(body, function(row)
          paste0("<tr>",
                 paste(sprintf("<td>%s</td>", row), collapse = ""),
                 "</tr>"),
          character(1)),
        collapse = "")
      style <- if (!is.null(max_width))
        sprintf(" style='max-width:%s;'", max_width) else ""
      HTML(paste0("<table class='coef-table'", style, ">",
                  header_html, body_html, "</table>"))
    }

    # Helper: signed LaTeX term, skipping zero coefficients.
    signed_term <- function(coef, expr) {
      if (abs(round(coef, 8)) < 1e-10) return("")
      if (coef >= 0) paste0(" + ", fmt(coef), expr)
      else            paste0(" - ", fmt(-coef), expr)
    }
    poly_latex <- function(i) {
      xterm <- paste0("(x - ", fmt(sp$x[i]), ")")
      paste0(fmt(sp$a[i]),
             signed_term(sp$b[i], xterm),
             signed_term(sp$c[i], paste0(xterm, "^2")),
             signed_term(sp$d[i], paste0(xterm, "^3")))
    }

    # ---- Step 1: sorted data points ----
    s1 <- tags$details(class = "step-box", open = NA,
      tags$summary(h4("Step 1: Sorted data points")),
      p(class = "step-friendly",
        "First, we line up your points from smallest to largest x."),
      mini_table(
        c("<em>i</em>", "<em>x<sub>i</sub></em>", "<em>y<sub>i</sub></em>"),
        lapply(seq_len(n + 1), function(i)
          c(as.character(i - 1), fmt(sp$x[i]), fmt(sp$y[i]))),
        max_width = "320px")
    )

    # ---- Step 2: interval widths ----
    s2 <- tags$details(class = "step-box", open = NA,
      tags$summary(h4(HTML("Step 2: Interval widths \\(h_i = x_{i+1} - x_i\\)"))),
      p(class = "step-friendly",
        "Next, we measure the gap between each pair of consecutive points."),
      mini_table(
        c("<em>i</em>",
          "<em>x<sub>i+1</sub></em> &minus; <em>x<sub>i</sub></em>",
          "<em>h<sub>i</sub></em>"),
        lapply(seq_len(n), function(i)
          c(as.character(i - 1),
            sprintf("%s &minus; %s", fmt(sp$x[i + 1]), fmt(sp$x[i])),
            fmt(sp$h[i]))),
        max_width = "380px")
    )

    # ---- Step 3: tridiagonal system ----
    if (n >= 2) {
      m <- nrow(sp$A)
      A_latex <- paste(
        apply(sp$A, 1, function(r) paste(sapply(r, fmt), collapse = " & ")),
        collapse = " \\\\ ")
      c_latex <- paste(sprintf("c_{%d}", seq_len(m)), collapse = " \\\\ ")
      r_latex <- paste(sapply(sp$rhs, fmt), collapse = " \\\\ ")
      matrix_eq <- sprintf(
        "$$\\begin{bmatrix} %s \\end{bmatrix}\\begin{bmatrix} %s \\end{bmatrix} = \\begin{bmatrix} %s \\end{bmatrix}$$",
        A_latex, c_latex, r_latex)
      s3 <- tags$details(class = "step-box", open = NA,
        tags$summary(h4(HTML("Step 3: Tridiagonal system \\(A\\,\\mathbf{c} = \\mathbf{r}\\) (interior \\(c_1 \\dots c_{n-1}\\))"))),
        p(class = "step-friendly",
          "To make the curve bend smoothly at every junction, ",
          "we set up a small system of equations."),
        p(HTML("$$h_i c_{i-1} + 2(h_i + h_{i+1})c_i + h_{i+1}c_{i+1} = 3\\left(\\frac{y_{i+1}-y_i}{h_{i+1}} - \\frac{y_i-y_{i-1}}{h_i}\\right)$$")),
        p(HTML("Natural boundary conditions: \\(c_0 = 0,\\ c_n = 0\\).")),
        p("Plugging your numbers into the formula above gives:"),
        p(HTML(matrix_eq))
      )
    } else {
      s3 <- tags$details(class = "step-box", open = NA,
        tags$summary(h4("Step 3: Tridiagonal system")),
        p("Only one interval — no interior c values to solve.")
      )
    }

    # ---- Step 4: solved c values ----
    s4 <- tags$details(class = "step-box", open = NA,
      tags$summary(h4(HTML("Step 4: Solve for \\(c_i\\)"))),
      p(class = "step-friendly",
        "Solving that system tells us how much each junction should curve."),
      mini_table(
        c("<em>i</em>", "<em>c<sub>i</sub></em>"),
        lapply(seq_len(n + 1), function(i)
          c(as.character(i - 1), fmt(sp$c_full[i]))),
        max_width = "260px")
    )

    # ---- Step 5: coefficients a, b, c, d ----
    s5 <- tags$details(class = "step-box solution-step", open = NA,
      tags$summary(h4(HTML("Step 5: Coefficients \\(a_i, b_i, c_i, d_i\\)"))),
      p(class = "step-friendly",
        "With those values in hand we can write down four numbers ",
        "(a, b, c, d) for every piece of the curve."),
      p(HTML("$$a_i = y_i,\\quad b_i = \\frac{y_{i+1}-y_i}{h_i} - \\frac{h_i(2c_i + c_{i+1})}{3},\\quad d_i = \\frac{c_{i+1}-c_i}{3 h_i}$$")),
      mini_table(
        c("<em>i</em>", "<em>x<sub>i</sub></em>",
          "<em>a<sub>i</sub></em>", "<em>b<sub>i</sub></em>",
          "<em>c<sub>i</sub></em>", "<em>d<sub>i</sub></em>"),
        lapply(seq_len(n), function(i)
          c(as.character(i - 1), fmt(sp$x[i]),
            fmt(sp$a[i]), fmt(sp$b[i]), fmt(sp$c[i]), fmt(sp$d[i]))))
    )

    # ---- Step 6: piecewise polynomials ----
    poly_blocks <- paste(
      sapply(seq_len(n), function(i)
        sprintf("$$S_{%d}(x) = %s, \\quad x \\in [%s, %s]$$",
                i - 1, poly_latex(i), fmt(sp$x[i]), fmt(sp$x[i + 1]))),
      collapse = "")
    s6 <- tags$details(class = "step-box solution-step", open = NA,
      tags$summary(h4("Step 6: Piecewise spline polynomials")),
      p(class = "step-friendly",
        "Each piece of the curve is now a tidy cubic polynomial — ",
        "here they are, one per interval."),
      HTML(poly_blocks)
    )

    # ---- Step 7: evaluation at x_query ----
    if (!is.na(xq) && xq >= sp$x[1] && xq <= sp$x[n + 1]) {
      i_r <- max(1, min(n, findInterval(xq, sp$x, all.inside = TRUE)))
      i   <- i_r - 1
      dx  <- xq - sp$x[i_r]
      val <- sp$a[i_r] + sp$b[i_r] * dx + sp$c[i_r] * dx^2 + sp$d[i_r] * dx^3
      s7 <- tags$details(class = "step-box", open = NA,
        tags$summary(h4(HTML(sprintf("Step 7: Evaluate at \\(x = %s\\)", fmt(xq))))),
        p(class = "step-friendly",
          "Finally, we plug your chosen x into the right piece ",
          "of the curve to read off the estimated y."),
        p(HTML(sprintf(
          "Since \\(%s \\in [%s,\\,%s]\\), we use piece \\(S_{%d}(x)\\).",
          fmt(xq), fmt(sp$x[i_r]), fmt(sp$x[i_r + 1]), i))),
        p(HTML(sprintf("$$\\Delta x = x - x_{%d} = %s - %s = %s$$",
                       i, fmt(xq), fmt(sp$x[i_r]), fmt(dx)))),
        p(HTML(sprintf(
          "$$S_{%d}(%s) = %s + (%s)(%s) + (%s)(%s)^2 + (%s)(%s)^3$$",
          i, fmt(xq),
          fmt(sp$a[i_r]),
          fmt(sp$b[i_r]), fmt(dx),
          fmt(sp$c[i_r]), fmt(dx),
          fmt(sp$d[i_r]), fmt(dx)))),
        div(class = "result-callout",
          HTML(sprintf("f(%s) = <span style='color:#F26D80;'>%s</span>",
                       fmt(xq), fmt(val))))
      )
    } else {
      s7 <- tags$details(class = "step-box", open = NA,
        tags$summary(h4("Step 7: Evaluation")),
        p(HTML(sprintf(
          "\\(x = %s\\) is outside \\([%s,\\,%s]\\) (extrapolation — not shown).",
          ifelse(is.na(xq), "?", fmt(xq)),
          fmt(sp$x[1]), fmt(sp$x[n + 1]))))
      )
    }

    withMathJax(tagList(s1, s2, s3, s4, s5, s6, s7))
  })

  # Output: presentation-friendly Solution view (MathJax)
  output$solution_ui <- renderUI({
    sp <- spline_data()
    xq <- input$x_query
    n  <- length(sp$h)

    fmt <- function(v) {
      r <- round(v, 5)
      if (length(r) == 1 && abs(r) < 1e-10) return("0")
      formatC(r, format = "g", digits = 5)
    }

    # Build one term of S_i(x) with proper sign and skip zero terms.
    signed_term <- function(coef, expr) {
      if (abs(round(coef, 8)) < 1e-10) return("")
      if (coef >= 0) paste0(" + ", fmt(coef), expr)
      else            paste0(" - ", fmt(-coef), expr)
    }

    poly_latex <- function(i) {
      xi <- sp$x[i]
      xterm <- paste0("(x - ", fmt(xi), ")")
      paste0(
        fmt(sp$a[i]),
        signed_term(sp$b[i], xterm),
        signed_term(sp$c[i], paste0(xterm, "^2")),
        signed_term(sp$d[i], paste0(xterm, "^3"))
      )
    }

    # --- Section 1: header ---
    header_html <- paste0(
      "<div class='card solution-section'>",
      "<h3 style='margin-top:0;'>Cubic Spline Solution</h3>",
      "<p style='color:#555;'>The natural cubic spline through your data, ",
      "shown as a coefficient table, a list of piecewise polynomials, ",
      "and a single combined piecewise function. Edit the data on the ",
      "<strong>Calculator</strong> tab and this view updates live.</p>",
      "</div>"
    )

    # --- Section 2: input data ---
    pts_rows <- paste(sapply(seq_len(n + 1), function(i) {
      sprintf("<tr><td>%d</td><td>%s</td><td>%s</td></tr>",
              i - 1, fmt(sp$x[i]), fmt(sp$y[i]))
    }), collapse = "")
    inputs_html <- paste0(
      "<div class='card solution-section'>",
      "<h4>Input Data</h4>",
      "<table class='coef-table' style='max-width:380px;'>",
      "<tr><th>i</th><th>x<sub>i</sub></th><th>y<sub>i</sub></th></tr>",
      pts_rows, "</table></div>"
    )

    # --- Section 3: coefficient table ---
    coef_rows <- paste(sapply(seq_len(n), function(i) {
      sprintf(
        "<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
        i - 1, fmt(sp$x[i]),
        fmt(sp$a[i]), fmt(sp$b[i]), fmt(sp$c[i]), fmt(sp$d[i])
      )
    }), collapse = "")
    coef_html <- paste0(
      "<div class='card solution-section'>",
      "<h4>Spline Coefficients</h4>",
      "<p>Each interval [x<sub>i</sub>, x<sub>i+1</sub>] has its own cubic",
      " $$S_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3.$$</p>",
      "<table class='coef-table'>",
      "<tr><th>i</th><th>x<sub>i</sub></th><th>a<sub>i</sub></th>",
      "<th>b<sub>i</sub></th><th>c<sub>i</sub></th><th>d<sub>i</sub></th></tr>",
      coef_rows, "</table></div>"
    )

    # --- Section 4: piecewise polynomial list (LaTeX) ---
    poly_lines <- sapply(seq_len(n), function(i) {
      sprintf("$$S_{%d}(x) = %s, \\quad x \\in [%s, %s]$$",
              i - 1, poly_latex(i), fmt(sp$x[i]), fmt(sp$x[i + 1]))
    })
    poly_html <- paste0(
      "<div class='card solution-section'>",
      "<h4>Piecewise Polynomials</h4>",
      paste(poly_lines, collapse = ""),
      "</div>"
    )

    # --- Section 5: combined piecewise function ---
    cases_rows <- paste(sapply(seq_len(n), function(i) {
      sprintf("%s, & x \\in [%s, %s]",
              poly_latex(i), fmt(sp$x[i]), fmt(sp$x[i + 1]))
    }), collapse = " \\\\ ")
    combined_html <- paste0(
      "<div class='card solution-section'>",
      "<h4>Combined Piecewise Function</h4>",
      "$$S(x) = \\begin{cases} ", cases_rows, " \\end{cases}$$",
      "</div>"
    )

    # --- Section 6: interpolated value callout ---
    if (is.na(xq)) {
      callout_html <- paste0(
        "<div class='solution-callout'>",
        "<strong>Enter a valid x_query on the Calculator tab.</strong></div>"
      )
    } else if (xq < sp$x[1] || xq > sp$x[n + 1]) {
      callout_html <- sprintf(paste0(
        "<div class='solution-callout' style='background:#fdecea; ",
        "border-color:#c44; border-left-color:#c44;'>",
        "<strong>Extrapolation:</strong> x_query = %s is outside the data ",
        "range [%s, %s]. The natural cubic spline is only defined within ",
        "the data range.</div>"),
        fmt(xq), fmt(sp$x[1]), fmt(sp$x[n + 1]))
    } else {
      i_r <- max(1, min(n, findInterval(xq, sp$x, all.inside = TRUE)))
      i   <- i_r - 1
      dx  <- xq - sp$x[i_r]
      val <- sp$a[i_r] + sp$b[i_r] * dx +
             sp$c[i_r] * dx^2 + sp$d[i_r] * dx^3
      callout_html <- sprintf(paste0(
        "<div class='solution-callout'>",
        "<div style='font-size:20px; margin-bottom:8px;'>",
        "<strong>f(%s) = %s</strong></div>",
        "Evaluated using <strong>S<sub>%d</sub></strong> because %s ",
        "&isin; [%s, %s].&nbsp; ",
        "&Delta;x = %s &minus; %s = %s.",
        "$$f(%s) = S_{%d}(%s) = %s$$",
        "</div>"),
        fmt(xq), fmt(val),
        i, fmt(xq), fmt(sp$x[i_r]), fmt(sp$x[i_r + 1]),
        fmt(xq), fmt(sp$x[i_r]), fmt(dx),
        fmt(xq), i, fmt(xq), fmt(val))
    }

    withMathJax(HTML(paste0(
      header_html, inputs_html, coef_html,
      poly_html, combined_html, callout_html
    )))
  })
}

# ---- Run App ----
shinyApp(ui, server)

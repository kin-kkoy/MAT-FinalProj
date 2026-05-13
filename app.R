# Cubic Spline Interpolation - R Shiny App
# Numerical Methods Final Project

library(shiny)
library(plotly)

# ---- Small UI helper: tooltip icon (uses native browser title attr) ----
help_tip <- function(text) {
  tags$span(class = "help-icon", title = text, HTML("&#9432;"))
}

# ---- Natural Cubic Spline (manual implementation) ----
# Solves the tridiagonal system for the second-derivative coefficients
# under natural boundary conditions (c_0 = c_n = 0), then derives
# a_i, b_i, c_i, d_i so each interval has
#   S_i(x) = a_i + b_i (x - x_i) + c_i (x - x_i)^2 + d_i (x - x_i)^3.
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

# ---- UI ----
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;600;700&display=swap")
  ),
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

    /* Title panel */
    .title-hero { padding: 28px 20px 6px; }
    .title-hero h1 { margin: 0; font-weight: 600; }
    .title-hero .tagline { color: #718096; margin-top: 6px; font-size: 16px; }

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
  ")),
  
  div(class = "title-hero",
      h1("Cubic Spline Interpolation"),
      div(class = "tagline",
          "Draw a smooth curve through any set of points — ",
          "an interactive guide to a classic numerical method.")
  ),
  
  tabsetPanel(id = "main_tabs",
              # ---- Introduction / Learn Tab ----
              tabPanel("Introduction",
                       fluidRow(
                         column(width = 10, offset = 1,
                                div(class = "card card-hero",
                                    h2("What is a cubic spline?"),
                                    div(class = "tag",
                                        "Imagine connecting dots on a graph — but instead of straight ",
                                        "lines that kink at every point, you want one perfectly smooth ",
                                        "curve. That's what a cubic spline does.",
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
                                    p(withMathJax("Between every pair of your data points, the spline lays down a small cubic polynomial — a curve of the form \\(a + bx + cx^2 + dx^3\\). Each piece is chosen so that where two pieces meet, they share the same value, the same slope, and the same curvature. That's why the joins are invisible.")),
                                    
                                    p(style = "color:#718096; font-size: 14px;",
                                      "The flavour used here is the ", tags$em("natural"), " cubic spline: at the two endpoints the curvature is set to zero, so the curve tapers off gently rather than flicking up."
                                    ),
                                    
                                    h3("Why does it work?"),
                                    
                                    p("Many real-world relationships are not linear. Cubic splines work because they enforce smooth connections between cubic segments, which naturally produces accurate approximations of real functions and minimizes errors that sharp or straight lines introduce. This makes later calculations (area, slope, motion) much more accurate and stable.")
                                    
                                    
                                ),
                                
                                div(class = "card",
                                    h3("Where is it used?"),
                                    tags$ul(class = "uses-list",
                                            tags$li(tags$strong("Animation & games — "),
                                                    "smooth motion paths for characters and cameras."),
                                            tags$li(tags$strong("Typography — "),
                                                    "the curves of letterforms in modern fonts."),
                                            tags$li(tags$strong("Computer-aided design — "),
                                                    "shaping car bodies, aircraft wings, product surfaces."),
                                            tags$li(tags$strong("Data analysis — "),
                                                    "filling in missing measurements between known samples."),
                                            tags$li(tags$strong("Engineering — "),
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
              ),
              
              # ---- Calculator Tab ----
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
              ),
              
              # ---- Solution Tab ----
              tabPanel("Solution",
                       fluidRow(
                         column(width = 10, offset = 1,
                                withMathJax(),
                                br(),
                                uiOutput("solution_ui")
                         )
                       )
              )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  default_x <- c(1, 2, 3, 4, 5, 6, 7)
  default_y <- c(2, 3, 5, 4, 6, 8, 7)
  n_points  <- reactiveVal(length(default_x))
  
  observeEvent(input$add_row,    n_points(n_points() + 1))
  observeEvent(input$remove_row, {
    if (n_points() > 3) n_points(n_points() - 1)
  })
  
  # Render the X/Y input table; preserve any values the user has typed in
  output$point_inputs <- renderUI({
    n <- n_points()
    isolate({
      get_val <- function(i, prefix, defaults) {
        v <- input[[paste0(prefix, "_", i)]]
        if (is.null(v)) defaults[((i - 1) %% length(defaults)) + 1] else v
      }
      cur_x <- vapply(seq_len(n), get_val, numeric(1),
                      prefix = "x", defaults = default_x)
      cur_y <- vapply(seq_len(n), get_val, numeric(1),
                      prefix = "y", defaults = default_y)
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
  
  # Reactive: gather and validate user-entered points
  data_points <- reactive({
    n <- n_points()
    read_col <- function(prefix) {
      vapply(seq_len(n), function(i) {
        v <- input[[paste0(prefix, "_", i)]]
        if (is.null(v)) NA_real_ else as.numeric(v)
      }, numeric(1))
    }
    x <- read_col("x")
    y <- read_col("y")
    
    validate(
      need(n >= 3, "Please provide at least 3 points."),
      need(all(!is.na(x)), "All X values must be filled in."),
      need(all(!is.na(y)), "All Y values must be filled in."),
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
    
    # Step 1 -- sorted data
    s1 <- div(class = "step-box",
              h4("Step 1: Sorted data points"),
              p(class = "step-friendly",
                "First, we line up your points from smallest to largest x."),
              tags$pre(paste(
                sprintf("(x_%d, y_%d) = (%s, %s)",
                        0:n, 0:n, fmt(sp$x), fmt(sp$y)),
                collapse = "\n"))
    )
    
    # Step 2 -- interval widths
    s2 <- div(class = "step-box",
              h4(withMathJax("Step 2: Interval widths $$h_i = x_{i+1} - x_i$$")),
              #h4("Step 2: Interval widths  h_i = x_{i+1} - x_i"),
              p(class = "step-friendly",
                "Next, we measure the gap between each pair of consecutive points."),
              tags$pre(paste(
                sprintf("h_%d = %s - %s = %s",
                        0:(n - 1), fmt(sp$x[2:(n + 1)]),
                        fmt(sp$x[1:n]), fmt(sp$h)),
                collapse = "\n"))
    )
    
    # Step 3 -- tridiagonal system
    if (n >= 2) {
      A_rows <- apply(sp$A, 1, function(r)
        paste(sprintf("%10s", fmt(r)), collapse = " "))
      rhs_row <- paste(sprintf("%10s", fmt(sp$rhs)), collapse = " ")
      s3 <- div(class = "step-box",
                h4("Step 3: Tridiagonal system  A · c = r  (interior c_1 ... c_{n-1})"),
                p(class = "step-friendly",
                  "To make the curve bend smoothly at every junction, ",
                  "we set up a small system of equations."),
                
                p(withMathJax("$$h_i c_{i-1} + 2(h_i + h_{i+1})c_i + h_{i+1}c_{i+1} = 3\\left(\\frac{y_{i+1}-y_i}{h_{i+1}} - \\frac{y_i-y_{i-1}}{h_i}\\right)$$")),
                p("Natural boundary conditions: c_0 = 0,  c_n = 0."),
                tags$pre(paste(c("A =", A_rows, "", "r =", rhs_row), collapse = "\n"))
      )
    } else {
      s3 <- div(class = "step-box",
                h4("Step 3: Tridiagonal system"),
                p("Only one interval — no interior c values to solve.")
      )
    }
    
    # Step 4 -- c values
    s4 <- div(class = "step-box",
              h4(withMathJax(paste0("Step 4: Solve for ", "\\(c_i\\)"))),
              p(class = "step-friendly",
                "Solving that system tells us how much each junction should curve."),
              tags$pre(paste(
                sprintf("c_%d = %s", 0:n, fmt(sp$c_full)),
                collapse = "\n"))
    )
    
    # Step 5 -- a, b, d
    rows <- paste(sprintf("  %d  %11s  %11s  %11s  %11s",
                          0:(n - 1), fmt(sp$a), fmt(sp$b),
                          fmt(sp$c), fmt(sp$d)), collapse = "\n")
    s5 <- div(class = "step-box solution-step",
              h4(withMathJax("Step 5: Coefficients \\(a_i, b_i, c_i, d_i\\)")),
              p(class = "step-friendly",
                "With those values in hand we can write down four numbers ",
                "(a, b, c, d) for every piece of the curve."),
              p(HTML("a<sub>i</sub> = y<sub>i</sub> &nbsp; b<sub>i</sub> = (y<sub>i+1</sub>-y<sub>i</sub>)/h<sub>i</sub> - h<sub>i</sub>(2c<sub>i</sub>+c<sub>i+1</sub>)/3 &nbsp; d<sub>i</sub> = (c<sub>i+1</sub>-c<sub>i</sub>)/(3h<sub>i</sub>)")),
              tags$pre(paste(
                "  i         a_i          b_i          c_i          d_i",
                rows, sep = "\n"))
    )
    
    # Step 6 -- spline polynomials
    poly_lines <- sapply(1:n, function(i) {
      sprintf("S_%d(x) = %s + %s (x - %s) + %s (x - %s)^2 + %s (x - %s)^3,    x in [%s, %s]",
              i - 1, fmt(sp$a[i]), fmt(sp$b[i]), fmt(sp$x[i]),
              fmt(sp$c[i]), fmt(sp$x[i]),
              fmt(sp$d[i]), fmt(sp$x[i]),
              fmt(sp$x[i]), fmt(sp$x[i + 1]))
    })
    s6 <- div(class = "step-box solution-step",
              h4("Step 6: Piecewise spline polynomials"),
              p(class = "step-friendly",
                "Each piece of the curve is now a tidy cubic polynomial — ",
                "here they are, one per interval."),
              tags$pre(paste(poly_lines, collapse = "\n"))
    )
    
    # Step 7 -- evaluation at x_query
    if (!is.na(xq) && xq >= sp$x[1] && xq <= sp$x[n + 1]) {
      i_r <- max(1, min(n, findInterval(xq, sp$x, all.inside = TRUE)))
      i   <- i_r - 1
      dx  <- xq - sp$x[i_r]
      val <- sp$a[i_r] + sp$b[i_r] * dx + sp$c[i_r] * dx^2 + sp$d[i_r] * dx^3
      s7 <- div(class = "step-box",
                h4(sprintf("Step 7: Evaluate at x_query = %s", fmt(xq))),
                p(class = "step-friendly",
                  "Finally, we plug your chosen x into the right piece ",
                  "of the curve to read off the estimated y."),
                tags$pre(paste(
                  sprintf("x_query = %s lies in [%s, %s], so use S_%d.",
                          fmt(xq), fmt(sp$x[i_r]), fmt(sp$x[i_r + 1]), i),
                  sprintf("Delta x = x_query - x_%d = %s - %s = %s",
                          i, fmt(xq), fmt(sp$x[i_r]), fmt(dx)),
                  sprintf("S_%d(%s) = %s + %s(%s) + %s(%s)^2 + %s(%s)^3",
                          i, fmt(xq),
                          fmt(sp$a[i_r]),
                          fmt(sp$b[i_r]), fmt(dx),
                          fmt(sp$c[i_r]), fmt(dx),
                          fmt(sp$d[i_r]), fmt(dx)),
                  sprintf("        = %s", fmt(val)),
                  sep = "\n"))
      )
    } else {
      s7 <- div(class = "step-box",
                h4("Step 7: Evaluation"),
                p(sprintf("x_query = %s is outside [%s, %s] (extrapolation — not shown).",
                          ifelse(is.na(xq), "NA", fmt(xq)),
                          fmt(sp$x[1]), fmt(sp$x[n + 1])))
      )
    }
    
    tagList(s1, s2, s3, s4, s5, s6, s7)
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

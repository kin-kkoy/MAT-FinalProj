# Cubic Spline Interpolation - R Shiny App
# Numerical Methods Final Project

library(shiny)

# ---- UI ----
ui <- fluidPage(
  titlePanel("Cubic Spline Interpolation"),

  tabsetPanel(
    # ---- Introduction Tab ----
    tabPanel("Introduction",
      fluidRow(
        column(width = 10, offset = 1,
          h3("What is Cubic Spline Interpolation?"),
          p("Cubic spline interpolation is a numerical method that fits a smooth",
            "curve through a set of data points by joining piecewise cubic",
            "polynomials between consecutive points."),
          p("Each cubic piece is constructed so that the overall curve is",
            "continuous, and its first and second derivatives are also continuous",
            "at every interior point. This makes the curve appear smooth without",
            "the oscillations often seen in high-degree polynomial interpolation."),
          h4("Why use it?"),
          tags$ul(
            tags$li("Produces a smooth curve through all given points."),
            tags$li("Avoids oscillations common with high-degree polynomials."),
            tags$li("Useful for estimating values between known data points.")
          ),
          h4("How to use this app"),
          tags$ol(
            tags$li("Go to the 'Calculator' tab."),
            tags$li("Enter your x-values and y-values as comma-separated numbers."),
            tags$li("Enter the x-value you want to interpolate."),
            tags$li("View the interpolated value and the spline plot.")
          )
        )
      )
    ),

    # ---- Calculator Tab ----
    tabPanel("Calculator",
      sidebarLayout(
        sidebarPanel(
          h4("Input Data Points"),
          textInput("x_vals", "X values (comma-separated):",
                    value = "1, 2, 3, 4, 5, 6, 7"),
          textInput("y_vals", "Y values (comma-separated):",
                    value = "2, 3, 5, 4, 6, 8, 7"),
          numericInput("x_query", "X to interpolate (x_query):", value = 3.5),
          helpText("Provide at least 3 points. X and Y must have the same length.")
        ),

        mainPanel(
          h4("Interpolated Value"),
          verbatimTextOutput("result"),
          h4("Spline Plot"),
          plotOutput("splinePlot", height = "450px")
        )
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # Parse comma-separated numeric input
  parse_nums <- function(text) {
    nums <- suppressWarnings(as.numeric(strsplit(text, ",")[[1]]))
    nums[!is.na(nums)]
  }

  # Reactive: parse and validate inputs
  data_points <- reactive({
    x <- parse_nums(input$x_vals)
    y <- parse_nums(input$y_vals)

    validate(
      need(length(x) >= 3, "Please provide at least 3 x-values."),
      need(length(y) >= 3, "Please provide at least 3 y-values."),
      need(length(x) == length(y), "X and Y must have the same number of values."),
      need(!any(duplicated(x)), "X values must be unique.")
    )

    # Sort by x to ensure proper interpolation
    ord <- order(x)
    list(x = x[ord], y = y[ord])
  })

  # Reactive: build the natural cubic spline function
  spline_fn <- reactive({
    d <- data_points()
    splinefun(d$x, d$y, method = "natural")
  })

  # Output: interpolated value at x_query
  output$result <- renderPrint({
    d <- data_points()
    f <- spline_fn()
    xq <- input$x_query

    validate(need(!is.na(xq), "Please enter a valid x_query."))

    yq <- f(xq)
    cat("f(", xq, ") = ", yq, "\n", sep = "")

    if (xq < min(d$x) || xq > max(d$x)) {
      cat("\nNote: x_query is outside the data range (extrapolation).")
    }
  })

  # Output: plot of points and spline curve
  output$splinePlot <- renderPlot({
    d <- data_points()
    f <- spline_fn()

    # Smooth curve over the data range
    x_seq <- seq(min(d$x), max(d$x), length.out = 500)
    y_seq <- f(x_seq)

    plot(d$x, d$y,
         pch = 19, col = "blue", cex = 1.3,
         xlab = "x", ylab = "y",
         main = "Cubic Spline Interpolation",
         xlim = range(c(d$x, input$x_query), na.rm = TRUE),
         ylim = range(c(d$y, y_seq), na.rm = TRUE))
    lines(x_seq, y_seq, col = "red", lwd = 2)

    # Mark the queried point if within plotting range
    xq <- input$x_query
    if (!is.na(xq)) {
      yq <- f(xq)
      points(xq, yq, pch = 4, col = "darkgreen", cex = 2, lwd = 3)
    }

    legend("topleft",
           legend = c("Data points", "Cubic spline", "Interpolated x_query"),
           col = c("blue", "red", "darkgreen"),
           pch = c(19, NA, 4),
           lty = c(NA, 1, NA),
           lwd = c(NA, 2, 3),
           bty = "n")
  })
}

# ---- Run App ----
shinyApp(ui, server)

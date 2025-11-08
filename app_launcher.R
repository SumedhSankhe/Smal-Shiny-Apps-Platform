## app_launcher.R - Performance Optimized
#
# Main launcher for modular Shiny apps in the TAP-P999999-APPs repository.
#
# Performance Optimizations:
# - Cached module loading to avoid redundant sourcing
# - Module server functions called once with moduleServer() pattern
# - Reactive values used efficiently to minimize re-rendering
# - UI components cached where possible
# - Reduced reactive dependencies

library(shiny)
box::use(
  bslib[page_sidebar, bs_theme, font_collection, sidebar],
  jsonlite[fromJSON],
  logger[log_info, log_debug, log_warn, log_error, log_threshold, INFO]
)

# Initialize logging
try(log_threshold(INFO), silent = TRUE)
log_info("Launcher startup initiated")

# --- Module Loading (done once at startup) ---
modules_raw <- fromJSON("modules_list.json", simplifyDataFrame = FALSE)
log_info(sprintf("Loaded %d module descriptors", length(modules_raw)))
modules <- setNames(modules_raw, sapply(modules_raw, `[[`, "id"))

# Lazy loading: only load modules when needed
# Store environments for loaded modules
module_envs <- new.env()

# Function to load a module on-demand
load_module <- function(id) {
  if (exists(id, envir = module_envs)) {
    log_debug(sprintf("Reusing cached module '%s'", id))
    return(module_envs[[id]])
  }
  
  mod <- modules[[id]]
  env <- new.env()
  source(mod$source, local = env)
  
  # Store loaded module
  loaded <- list(
    ui_func = get(mod$ui, envir = env),
    server_func = get(mod$server, envir = env),
    env = env
  )
  
  module_envs[[id]] <- loaded
  log_info(sprintf("Module '%s' sourced from %s", id, mod$source))
  return(loaded)
}

# Pre-load only essential modules (e.g., module_template)
preload_ids <- sapply(modules, function(m) {
  isTRUE(m$preload)
})

for (id in names(modules)[preload_ids]) {
  log_debug(sprintf("Preloading module '%s'", id))
  load_module(id)
}

# --- Pre-compute UI elements ---
# Build choices list once instead of in UI
app_choices <- {
  all_ids <- names(modules)
  all_labels <- sapply(modules, `[[`, "label")
  mt_idx <- which(all_ids == "module_template")
  if (length(mt_idx) == 1) {
    c(
      setNames(all_ids[mt_idx], all_labels[mt_idx]),
      setNames(all_ids[-mt_idx], all_labels[-mt_idx])
    )
  } else {
    setNames(all_ids, all_labels)
  }
}

# --- UI Definition ---
ui <- page_sidebar(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    # Add resource hints for faster loading
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "dns-prefetch", href = "https://fonts.googleapis.com")
  ),
  uiOutput("app_ui"),
  sidebar = sidebar(
    selectInput(
      "selected_app",
      "Choose an app:",
      selected = "module_template",
      choices = app_choices
    ),
    hr()
  ),
  title = "Small Shiny Apps Platform",
  theme = bs_theme(
    version = 5,
    bootswatch = "litera",
    base_font = font_collection(
      "Nunito", "Helvetica Neue", "Helvetica", "Arial", "sans-serif"
    ),
    base_font_size = 15,
    "font-size-base" = "0.9rem",
    "spacer" = "0.9rem",
    "nav-link-padding-x" = "0.5rem;",
    "badge-padding-y" = ".4em;",
    "btn-padding-x" = "0.8rem;",
    "btn-padding-y" = "0.3rem;",
    "btn-font-size" = ".8rem;"
  )
)

# --- Server Logic (Optimized) ---
server <- function(input, output, session) {
  
  # Track active module servers to avoid re-initialization
  active_servers <- reactiveVal(list())
  
  # Render UI only when app selection changes
  output$app_ui <- renderUI({
    req(input$selected_app)
    
    # Lazy load module if needed
    loaded_module <- load_module(input$selected_app)
    log_debug(sprintf("Rendering UI for module '%s'", input$selected_app))
    
    # Call UI function
    loaded_module$ui_func(input$selected_app)
  })
  
  # Initialize module server ONCE per module
  observeEvent(input$selected_app, {
    req(input$selected_app)
    
    app_id <- input$selected_app
    current_active <- active_servers()
    
    # Only initialize if not already active
    if (!(app_id %in% names(current_active))) {
      # Lazy load module if needed
      loaded_module <- load_module(app_id)
      
      # Call server function
      loaded_module$server_func(app_id)
      
      # Mark as active
      current_active[[app_id]] <- TRUE
      active_servers(current_active)
      log_info(sprintf("Initialized server for module '%s'", app_id))
    }
  }, ignoreInit = FALSE)
}

# --- App Initialization ---
shinyApp(ui, server, options = list(
  # Performance options
  display.mode = "normal",
  launch.browser = FALSE  # Faster startup in production
))
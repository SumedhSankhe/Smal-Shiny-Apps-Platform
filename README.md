# TAP-Utility-WebApp

Welcome! This is a collection of small, easy-to-use apps for lab and data work. You can run any app by itself, or use the main launcher to pick from a menu.

> **Note:** This application is now publicly accessible. You can run it locally via `app_launcher.R` or access a hosted deployment if provided.

------------------------------------------------------------------------

## 🚀 Quick Start

1.  **Open RStudio**.
2.  **Open the file** `app_launcher.R`.
3.  Click the **Run App** button at the top of RStudio.
4.  In the app window, **choose a tool** from the dropdown menu.
5.  Use the app! If you get stuck, see the FAQ below or contact support.

------------------------------------------------------------------------

## ❓ What Can I Do With This?

-   Quickly make plots or randomize samples for lab work
-   Use simple tools for common data tasks
-   No coding needed—just click and go!

------------------------------------------------------------------------

## ⚡ FAQ (Frequently Asked Questions)

**Q: I see an error or the app won’t start!** A: Close RStudio and try again. If it still doesn’t work, contact support (see below).

**Q: Can I use two apps together or share data between them?** A: No, each app is separate for safety. If you need this, contact support.

**Q: I want a new tool or feature!** A: Great! Email your idea to support.

------------------------------------------------------------------------

## 🆘 Who to Contact for Help

If you have any questions, problems, or want a new feature, please contact:

**Sumedh Sankhe**\
Email: [sumedh.sankhe@gmail.com](mailto:sumedh.sankhe@gmail.com)

------------------------------------------------------------------------

## 📚 Glossary

-   **App**: A small tool you can run from the menu.
-   **Launcher**: The main menu that lets you pick which app to use.
-   **Module**: (For developers) A technical term for each app/tool.
-   **RStudio**: The program you use to run these apps.

------------------------------------------------------------------------

## For Developers & Advanced Users

### Folder Structure & Module Organization

``` text
TAP-Utility-WebApp/
├── app_launcher.R           # Main launcher app
├── modules_list.json        # Registry of available modules
├── R/                      # All modules live here, each in its own folder
│   ├── Module_Template/    # Template for new modules
│   │   └── app.R
│   └── New_Module/         # (Example) Your new module folder
│       └── app.R
├── renv/                   # Project-local R package library
├── README.md
└── ...
```

### How to Add a New Module

1.  Copy `R/Module_Template` to `R/Your_New_Module`.
2.  Rename the UI/server functions in `app.R` to match your module.
3.  Add any data or helper scripts inside your module folder.
4.  Register your module in `modules_list.json` (see below).

**Best Practices:** - Keep each module self-contained for easy reuse and maintenance. - Use the `{box}` package for explicit imports (see below). - Document your module at the top of its `app.R` file.

### Caution: Data File Size Limits

**Do NOT include large data files or any files totaling more than 50MB in this repository.**

-   Large files can cause problems with version control, sharing, and deployment.
-   If your app requires large data, use a small sample or synthetic dataset for development, and document how to obtain or link to the real data externally.
-   Use processed and sanitized data stored as .rds or .rda files

### Using the {box} Package in Modules

All modules in this repository should use the [{box}](https://klmr.me/box/) package for explicit, modular imports instead of `library()` or `::` calls. This ensures that dependencies are clear, only the needed functions are imported, and there is no global namespace pollution.

**Note:** The `shiny` package is provided in the global scope by the launcher platform, so you do not strictly need to import `shiny` with `box::use` in your module. However, using `box::use` for `shiny` is recommended for explicitness and consistency, especially if you want to test your module standalone or reuse it elsewhere.

**How to use {box} in your module:**

1.  At the top of your `app.R` (or module file), import only the functions you need from each package:

``` r
box::use(
  shiny[fluidPage, NS, textInput, textOutput, renderText, ...],
  bslib[layout_column_wrap, card, card_title, card_body],
  markdown[markdownToHTML]
)
```

2.  When calling imported functions inside your module, use them directly (e.g., `fluidPage(...)`), **except** inside a function body (like `server`) where you must use the namespace (e.g., `shiny$renderText(...)`).

-   This is required because of how {box} handles scoping.

3.  Do **not** use `library()` or `::` in your module code.
4.  See `R/Module_Template/app.R` for a complete example.

**Best practices:** - Only import the functions you actually use. - Use the `shiny$` prefix for Shiny functions inside server functions. - Keep all imports at the top of your file for clarity.

For more details, see the [box documentation](https://klmr.me/box/).

### Editing the Module List

-   Open `modules_list.json` in any text editor.

-   Each module entry should look like:

    ``` json
    {
      "id": "your_module_id",
      "label": "Your Module Name",
      "source": "R/Your_Module_Folder/app.R",
      "ui": "your_module_ui",
      "server": "your_module_server"
    }
    ```

-   The `id` must be unique and match the folder/function naming.

### Requirements

-   R (\>= 4.0)
-   R packages: shiny, jsonlite, (and any packages required by your modules)
-   Use `renv` to manage dependencies: `renv::restore()`
-   For using any additional packages that are required use the `renv::install("packageName", lock = TRUE)`

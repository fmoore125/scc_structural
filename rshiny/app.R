## setwd("~/research/scciams/scc_structural/rshiny")

source("../src/analysis/randomforest_dists_lib.R", chdir=T)
load("../outputs/rfdistsmodel-final-precompute.RData")
load("dats.RData")

library(shiny)
library(shinyjs)
library(plotly)
library(tibble)

strchoices <- list("Ambiguity/Model Uncertainty"="Ambiguity/Model Uncertainty",
                   "Earth System"="Earth system",
                   "Epstein-Zin"="Epstein-Zin",
                   "Inequality Aversion"="Inequality Aversion",
                   "Learning"="Learning",
                   "Limited Substitutability"="Limitedly-Substitutable Goods",
                   "Persistent/Growth Damages"="Persistent / Growth Damages",
                   "Tipping Points: Climate"="Tipping Points",
                   "Tipping Points: Damages"="Tipping Points2")
strchoices2 <- list()
for (nn in names(strchoices)) {
    strchoices2[[paste0(nn, " (", round(100 * mean(dat[, strchoices[[nn]]] == "1")), "% of lit.; ",
                        round(100 * mean(idealdat[, strchoices[[nn]]] == "1")), "% by exp.)")]] <-
        strchoices[[nn]]
}

uncchoices <- list()
for (cc in which(names(dat) == 'TFP Growth'):which(names(dat) == 'Risk Aversion (EZ Utility)'))
    uncchoices[[paste0(names(dat)[cc], " (", round(100 * mean(dat[, cc] == "1")), "% of lit.)")]] <- names(dat)[cc]

othchoices <- c("Backstop Price?", "Declining Discounting?",
                "Market Only Damages", "Other Market Failure?")
othchoices2 <- list()
for (nn in othchoices) {
    othchoices2[[paste0(nn, " (", round(100 * mean(dat[, nn] == "1")), "% of lit.)")]] <- nn
}

ui <- fluidPage(
    useShinyjs(),
    tags$head(
             tags$style(HTML("
            /* Custom CSS for sidebar scrolling */
            .sidebar {
                position: fixed;
                height: 100%;
                overflow-y: auto;
            }
            .config-table {
                width: 100%;
                border-collapse: collapse;
            }
            .config-table th, .config-table td {
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }
            .config-table th {
                background-color: #f2f2f2;
            }
            .config-table tr:nth-child(even) {
                background-color: #f9f9f9;
            }
        "))
        ),
    titlePanel("Modeling the SCC"),
    sidebarLayout(
        sidebarPanel(
            class='sidebar',
            plotlyOutput("discountdist", height="200px"),
            sliderInput("selected_rate",
                        "Select Discount Rate:",
                        min=min(dat$discountrate, na.rm=T),
                        max=round(max(dat$discountrate, na.rm=T)),
                        value=mean(idealdat$discountrate, na.rm=T),
                        step=0.1),
            plotlyOutput("damagedist", height="200px"),
            sliderInput("selected_dscc",
                        "Select Damage-based SCC:",
                        min=0,
                        max=round(quantile(idealdat$log.scc.synth, .99, na.rm=T)),
                        value=mean(idealdat$log.scc.synth, na.rm=T),
                        step=0.1),
            checkboxGroupInput(
                inputId="structgroup",
                label="Choose structural changes:",
                choices=strchoices2),
            checkboxGroupInput(
                inputId="uncgroup",
                label="Choose uncertainty:",
                choices=uncchoices),
            checkboxGroupInput(
                inputId="othergroup",
                label="Choose other options:",
                choices=othchoices2)
        ),
        mainPanel(
            plotlyOutput("sccdist", height="500px"),
            uiOutput("configTable"),
            hidden(actionButton("saveButton", "Save Configuration")),
            hidden(actionButton("loadButton", "Load Selected Configuration"))
        )
    )
)

makeConfigRow <- function(rate, dscc, quants) {
    midpoints <- (quants[-length(quants)] + quants[-1]) / 2
    mu <- sum(diff(all.qs) * midpoints)

    tibble(`Discount Rate`=paste0(rate, "%"), `Damage-based SCC`=paste0("exp(", dscc, ")"),
           `Mean SCC`=round(mu), `Median SCC`=round(quants[all.qs == .5]),
           `IQR`=paste(round(quants[all.qs == .25]), "-", round(quants[all.qs == .75])))
}

constructTableHTML <- function(rows) {
    if (!is.null(rows)) {
        configRows <- apply(rows, 1, function(row) {
            paste0("<tr>",
                   paste(sapply(row, function(cell) {
                       paste0("<td>", cell, "</td>")
                   }), collapse=""),
                   "</tr>")
        })
        tableHeader <- paste0("<tr>",
                              paste(sapply(colnames(rows), function(col) {
                                  paste0("<th>", col, "</th>")
                              }), collapse=""),
                              "</tr>")
        tableHTML <- paste0("<table class='config-table'>",
                            tableHeader,
                            paste(configRows, collapse=""),
                            "</table>")

        HTML(tableHTML)
    }
}

server <- function(input, output, session) {
    output$discountdist <- renderPlotly({
        plot_ly(height=200) %>%
            add_histogram(x=dat$discountrate, opacity=.6, name="Literature Distribution") %>%
            add_histogram(x=idealdat$discountrate, opacity=.6, name="Drupp et al. Distribution") %>%
            add_trace(x=4.6, y=0, type="scatter", mode="markers", marker=list(color="blue", size=10), name="DICE discount rate") %>%
            layout(title="Distribution of Discount Rates",
                   xaxis=list(title="Discount Rate (%)"),
                   yaxis=list(title=NA), #yaxis=list(title="Relative Frequency"),
                   legend=list(x=1, y=1, xanchor='right', yanchor='top'),
                   barmode='overlay', margin=list(l=30, r=0, b=0)) %>%
            add_trace(x=c(input$selected_rate), type="scatter", mode="markers",
                      marker=list(color="red", size=10),
                      name="Selected Rate")
    })
    output$damagedist <- renderPlotly({
        allowed.lss <-
        plot_ly(height=200) %>%
            add_histogram(x=idealdat$log.scc.synth[idealdat$log.scc.synth >= 0 & idealdat$log.scc.synth <= quantile(idealdat$log.scc.synth, .99, na.rm=T)], opacity=.6, name="Literature Distribution") %>%
            add_histogram(x=dicedat$log.scc.synth[dicedat$log.scc.synth >= 0 & dicedat$log.scc.synth <= quantile(idealdat$log.scc.synth, .99, na.rm=T)], opacity=.6, name="DICE Distribution") %>%
            layout(title="Distribution of Damage-based SCCs",
                   xaxis=list(title="Log of Damage-based SCC (log($/t))"),
                   yaxis=list(title=NA),
                   legend=list(x=1, y=1, xanchor='right', yanchor='top'),
                   barmode='overlay', margin=list(l=30, r=0, b=0)) %>%
            add_trace(x=c(input$selected_dscc), type="scatter", mode="markers",
                      marker=list(color="red", size=10),
                      name="Selected DSCC")
    })

    calcSCCValues <- reactive({
        preddat <- data.frame(discountrate=input$selected_rate, log.scc.synth=input$selected_dscc,
                              sccyearformerge=2020, Year=mean(dat$Year))
        for (nn in names(strchoices))
            preddat[, strchoices[[nn]]] <- ifelse(strchoices[[nn]] %in% input$structgroup, "1", "0")
        for (nn in names(uncchoices))
            preddat[, uncchoices[[nn]]] <- ifelse(uncchoices[[nn]] %in% input$uncgroup, "1", "0")
        for (nn in othchoices)
            preddat[, nn] <- ifelse(nn %in% input$othergroup, "1", "0")
        predict.forest(forest, preddat, NULL)
    })

    output$sccdist <- renderPlotly({
        output <- calcSCCValues()
        xx <- c()
        yy <- c()

        for (ii in 2:length(all.qs)) {
            xx <- c(xx, output[ii-1], output[ii])
            yy <- c(yy, rep((all.qs[ii] - all.qs[ii-1]) / (output[ii] - output[ii-1]), 2))
        }
        xx <- c(xx, output[length(output)], output[1])
        yy <- c(yy, 0, 0)

        pp <- plot_ly(height=500)

        if (!is.null(savedConfigs$outputs)) {
            xx2 <- c()
            yy2 <- c()
            group2 <- c()

            for (jj in 1:length(savedConfigs$outputs)) {
                output <- savedConfigs$outputs[[jj]]
                for (ii in 2:length(all.qs)) {
                    xx2 <- c(xx2, output[ii-1], output[ii])
                    yy2 <- c(yy2, rep((all.qs[ii] - all.qs[ii-1]) / (output[ii] - output[ii-1]), 2))
                }
                xx2 <- c(xx2, output[length(output)], output[1])
                yy2 <- c(yy2, 0, 0)
                group2 <- c(group2, rep(paste('Saved', jj), length(all.qs) * 2))
            }

            pp <- pp %>%
                add_polygons(x=xx2, y=yy2, split=group2, fillcolor='rgba(0,0,0,0)', line=list(width=1, color='toself'), name=group2)
        }

        pp %>%
            add_polygons(x=xx, y=yy, split="Current", fillcolor="blue", line=list(width=1, color='blue'), name="Current") %>%
            layout(title="Distribution of Social Cost of Carbon",
                   xaxis=list(title="SCC ($/t)"),
                   yaxis=list(title="Probability Density"),
                   legend=list(x=1, y=1, xanchor='right', yanchor='top'))
    })

    ## Define a reactiveValues object to store saved configurations
    savedConfigs <- reactiveValues(data=NULL, rows=NULL, outputs=NULL)

    ## Table to display saved configurations
    output$configTable <- renderUI({
        output <- calcSCCValues()

        constructTableHTML(rbind(cbind(makeConfigRow(input$selected_rate, input$selected_dscc, output),
                                       action=paste0('<button class="btn btn-success" onclick="$(\'#saveButton\').click(); return false;">Save</button>')),
                                 savedConfigs$rows))
    })

    ## Save configuration on button click
    observeEvent(input$saveButton, {
        newEntry <- data.frame(
            discountrate=input$selected_rate,
            dscc=input$selected_dscc,
            structural=paste(input$structgroup, collapse=", "),
            uncertainty=paste(input$uncgroup, collapse=", "),
            others=paste(input$othergroup, collapse=", ")
        )

        output <- calcSCCValues()
        newRow <- cbind(makeConfigRow(input$selected_rate, input$selected_dscc, output),
                        action=paste0('<button class="btn btn-success" onclick="Shiny.onInputChange(\'loadButton\', this.parentElement.parentElement.rowIndex - 1);">Load</button>'))

        if (is.null(savedConfigs$data)) {
            savedConfigs$data <- newEntry
            savedConfigs$rows <- newRow
            savedConfigs$outputs <- list(output)
        } else {
            savedConfigs$data <- rbind(savedConfigs$data, newEntry)
            savedConfigs$rows <- rbind(savedConfigs$rows, newRow)
            savedConfigs$outputs[[length(savedConfigs$outputs)+1]] <- output
        }
    })

    ## Load configuration on button click
    observeEvent(input$loadButton, {
        selectedRow <- as.numeric(input$loadButton)
        if (!is.null(savedConfigs$data) && selectedRow <= nrow(savedConfigs$data)) {
            updateSliderInput(session, "selected_rate", value=savedConfigs$data[selectedRow, "discountrate"])
            updateSliderInput(session, "selected_dscc", value=savedConfigs$data[selectedRow, "dscc"])
            updateCheckboxGroupInput(session, "structgroup", selected=strsplit(savedConfigs$data[selectedRow, "structural"], ", ")[[1]])
            updateCheckboxGroupInput(session, "uncgroup", selected=strsplit(savedConfigs$data[selectedRow, "uncertainty"], ", ")[[1]])
            updateCheckboxGroupInput(session, "othergroup", selected=strsplit(savedConfigs$data[selectedRow, "others"], ", ")[[1]])
        }
    })
}

shinyApp(ui=ui, server=server)

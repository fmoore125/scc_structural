## setwd("~/research/scciams/scc_structural/rshiny")

source("randomforest_dists_lib.R", chdir=T)
load("rfdistsmodel-final-precompute.RData")
load("dats.RData")

library(shiny)
library(shinyjs)
library(plotly)
library(tibble)

addResourcePath(prefix='static', directoryPath='static')

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
    strchoices2[[paste0(nn, " (", signif(100 * mean(dat[, strchoices[[nn]]] == "1"), 2), "% of lit.; ",
                        signif(100 * mean(idealdat[, strchoices[[nn]]] == "1"), 2), "% by exp.)")]] <-
        strchoices[[nn]]
}

uncchoicechange <- list("Carbon Cycle2"="Carbon Cycle",
                        "EMUC2"="Elasticity of marginal utility",
                        "PRTP2"="Pure rate of time preference",
                        "Risk Aversion (EZ Utility)"="Risk Aversion")

uncchoices <- list()
for (cc in which(names(dat) == 'TFP Growth'):which(names(dat) == 'Risk Aversion (EZ Utility)')) {
    label <- names(dat)[cc]
    if (label %in% names(uncchoicechange))
        label <- uncchoicechange[[label]]
    uncchoices[[paste0(label, " (", signif(100 * mean(dat[, cc] == "1"), 2), "% of lit.)")]] <- names(dat)[cc]
}

othchoices <- c("Backstop Price?", "Declining Discounting?",
                "Market Only Damages", "Other Market Failure?")
othchoices2 <- list()
for (nn in othchoices) {
    othchoices2[[paste0(nn, " (", signif(100 * mean(dat[, nn] == "1"), 2), "% of lit.)")]] <- nn
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
    tags$script(HTML('
const tooltips = {
    // Structural changes
    "Ambiguity/Model Uncertainty": "Considers various models showing pessimistic results for SCC, leading to cautionary decision-making and higher SCC values.",
    "Earth system": "Changes to how emissions translate into atmospheric concentrations and temperatures, impacting SCC based on temperature response.",
    "Epstein-Zin": "Incorporates separated risks across time and nature, influencing SCC by preventing higher risk aversion from lowering discount rates.",
    "Inequality Aversion": "Weighting impacts to reflect higher welfare losses in poorer regions affected more by climate change, often increasing the SCC.",
    "Learning": "Models changing over time as data updates unknown parameters, usually reducing SCC by refining policy to match true climate states.",
    "Limitedly-Substitutable Goods": "Models considering environmental goods and their substitutability with consumption goods, with ambiguous effect on SCC.",
    "Persistent / Growth Damages": "Accounts for long-lasting effects of temperature changes on productivity or capital depreciation, leading to higher SCC.",
    "Tipping Points": "Incorporates subsystems like ice sheets or ocean circulation that can suddenly change state, affecting the SCC diversely.",
    "Tipping Points2": "Models sudden, irreversible economic impacts from climate thresholds, generally increasing SCC due to higher future damages.",

    // Parametric uncertainty
    "TFP Growth": "Uncertainty in Total Factor Productivity growth rate.",
    "Population Growth": "Variability in future population growth predictions.",
    "Emissions Growth": "Uncertainty in future emissions growth rates.",
    "Transient Climate Response": "Short-term temperature sensitivity to CO2 emissions.",
    "Carbon Cycle2": "Uncertainty in carbon cycle dynamics affecting atmospheric CO2.",
    "Equilibrium Climate Sensitivity": "Long-term temperature sensitivity to CO2 levels.",
    "Tipping Point Magnitude": "Expected impact size of crossing climate or damage tipping points.",
    "Damage Function": "Uncertainty in the relationship between climate changes and economic damages.",
    "Adaptation Rates": "Rates at which populations or economies adapt to climate impacts.",
    "Income Elasticity": "Relation between income changes and changes in climate damage resilience or vulnerability.",
    "Constant Discount Rate": "Uncertainty in the fixed rate used to discount future costs and benefits to present value.",
    "EMUC2": "Elasticity of marginal utility of consumption, relating to diminishing returns on consumption increases.",
    "PRTP2": "Pure rate of time preference, reflecting how future utility is discounted.",
    "Risk Aversion (EZ Utility)": "Degree to which variability in outcomes affects welfare valuation.",

    // Other options
    "Backstop Price?": "Consideration of a future technology eliminating emissions at a fixed cost.",
    "Declining Discounting?": "Uses discount rates that decrease over time instead of remaining constant.",
    "Market Only Damages": "Considers only direct market damages, excluding non-market and indirect effects.",
    "Other Market Failure?": "Accounts for additional market inefficiencies beyond carbon pricing."
};

  $(document).ready(function() {
    $("input:checkbox").each(function() {
      $(this).prop({checked: false, indeterminate: true});
      var key = $(this).attr("value");
      var $sister = $(":checkbox[name=" + $(this).prop("name") + "_deter][value=\'" + key + "\']");
      $sister.prop("checked", false);
      $(this).on("change", function() {
        if ($(this).prop("checked") && !$sister.prop("checked")) {
          // Transition from minus to plus
          $sister.prop("checked", true);
          $sister.trigger("change");
          return true;
        } else if ($(this).prop("checked") && $sister.prop("checked")) {
          // Transition from blank to minus
          $(this).prop("checked", false);
          $(this).prop("indeterminate", true);
          $sister.prop("checked", false);
          $sister.trigger("change");
          return true;
        } else {
          // Transition from plus to blank
          $sister.prop("checked", true);
          $sister.trigger("change");
          return true;
        }
      });

      var tooltip = $("<img style=\\"margin-left: 10px\\" src=\\"static/dialog-question-symbolic.svg\\"/>");
      $(tooltip).appendTo($(this).parent()).attr("title", tooltips[key]).tooltip();
    });
  });
')),
    titlePanel("Social Cost of Carbon Generator"),
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
            hidden(
                checkboxGroupInput(
                    inputId="structgroup_deter",
                    label="Choose structural changes:",
                    choices=strchoices2)),
            checkboxGroupInput(
                inputId="uncgroup",
                label="Choose uncertainty:",
                choices=uncchoices),
            hidden(
                checkboxGroupInput(
                    inputId="uncgroup_deter",
                    label="Choose uncertainty:",
                    choices=uncchoices)),
            checkboxGroupInput(
                inputId="othergroup",
                label="Choose other options:",
                choices=othchoices2),
            hidden(
                checkboxGroupInput(
                    inputId="othergroup_deter",
                    label="Choose other options:",
                    choices=othchoices2)),
            HTML('<div style="height: 40px"></div>')
        ),
        mainPanel(
            plotlyOutput("sccdist", height="500px"),
            uiOutput("configTable"),
            hidden(actionButton("saveButton", "Save Configuration")),
            hidden(actionButton("loadButton", "Load Selected Configuration")),
            HTML("Note: This tool shows predictions from a random forest model trained on published SCC estimates. SCC distributions from different model structures and parameter values are based on subsets of the literature and should be interpreted with caution.")
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

get.preddat <- function(input) {
    print(input)
    preddat <- data.frame(discountrate=input$selected_rate, log.scc.synth=input$selected_dscc,
                          sccyearformerge=2020, Year=mean(dat$Year))
    ## Only add if the checkbox is clicked
    print(input$structgroup_deter)
    print(input$structgroup)
    for (nn in names(strchoices))
        if (strchoices[[nn]] %in% input$structgroup_deter)
            preddat[, strchoices[[nn]]] <- ifelse(strchoices[[nn]] %in% input$structgroup, "1", "0")
    for (nn in names(uncchoices))
        if (uncchoices[[nn]] %in% input$uncgroup_deter)
            preddat[, uncchoices[[nn]]] <- ifelse(uncchoices[[nn]] %in% input$uncgroup, "1", "0")
    for (nn in othchoices)
        if (nn %in% input$othergroup_deter)
            preddat[, nn] <- ifelse(nn %in% input$othergroup, "1", "0")
    preddat
}

## preddat = get.preddat(list(selected_rate=2, selected_dscc=2))

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

    saved.preddat <- reactive({
        get.preddat(input)
    })

    calcSCCValues <- reactive({
        preddat <- saved.preddat()

        predict.forest(forest, preddat, NULL)
    })

    output$sccdist <- renderPlotly({
        preddat <- saved.preddat()
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

library(shiny)
library(shinybrowser)
library(rhandsontable)

library(ggplot2)
library(plotly)

# Define server logic
shinyServer(function(input, output) {
    
    # Normalize utility cost associated with a false positive
    d_01 = 1
    # Range of utility costs associated with a false negative
    delta_10_vect = 10^seq(-1, 1, by = 0.1)
    # Test utility
    K = 0
    
    # Initialize an example table of values for a 5-gene breast panel
    table_of_values = data.frame(
        Gene = c("ATM", "BRCA1", "BRCA2", "CHEK2", "PALB2"), 
        Prevalence = c(0.00380, 0.00117, 0.00135, 0.00519, 0.00114), 
        Penetrance = c(0.34513, 0.73206, 0.71665, 0.19461, 0.38481),
        Precision = c(1e+02, 1e+04, 1e+04, 1e+02, 1e+02),
        stringsAsFactors = FALSE) 
    
    # Store reactive values
    values = reactiveValues()
    
    # Convert initial handsontable to R data frame
    observe({
        if (!is.null(input$table_of_values)) {
            table_of_values = hot_to_r(input$table_of_values)
        } else {
            if (is.null(values[["table_of_values"]]))
                table_of_values <- table_of_values
            else
                table_of_values <- values[["table_of_values"]]
        }
        values[["table_of_values"]] <- table_of_values
    })
    
    # Format table appearance
    output$table_of_values <- renderRHandsontable({
        table_of_values <- values[["table_of_values"]]
        if (!is.null(table_of_values))
            rhandsontable(table_of_values, useTypes = TRUE, 
                          stretchH = "all", maxRows = 10) %>% # At most 10 rows
            hot_col(c("Prevalence", "Penetrance"), format = "0.00000") %>% # Round to 5 decimals
            hot_col("Precision", format = "0") %>% # Round to closest integer
            hot_context_menu(allowRowEdit = TRUE) %>% # Allow user to add/remove rows
            hot_cols(columnSorting = TRUE) # Allow user to sort columns
    })
    
    # Pull in handsontable (possibly updated by user) to a reactive R data frame
    table_of_values2 <- reactive({
        if (!is.null(input$table_of_values)) {
            hot_to_r(input$table_of_values)
        } else {
            if (is.null(values[["table_of_values"]]))
                table_of_values
            else
                values[["table_of_values"]]
        }
    })
    
    # Generate plots based on user inputs
    output$plots <- renderPlotly({
        # Table of values
        table_of_values2 = table_of_values2()
        # Drop all empty rows and rows with missing values
        table_of_values2[is.null(table_of_values2)] = NA
        table_of_values2 = na.omit(table_of_values2)
        
        # Save columns of table as vectors
        genes = as.character(table_of_values2$Gene)
        prev = table_of_values2$Prevalence
        pen = table_of_values2$Penetrance
        precision = table_of_values2$Precision
        
        # Check that prevalences, penetrances, and precisions are within valid range
        if (any(prev < 0) || any(prev > 1)) {
            stop(paste("Prevalences out of bounds for the following genes:", 
                       paste(genes[prev < 0 | prev > 1], collapse = ", ")))
        }
        if (any(pen < 0) || any(pen > 1)) {
            stop(paste("Penetrances out of bounds for the following genes:", 
                       paste(genes[pen < 0 | pen > 1], collapse = ", ")))
        }
        if (any(precision <= 0)) {
            stop(paste("Non-positive precisions for the following genes:", 
                       paste(genes[precision <= 0], collapse = ", ")))
        }
        
        # Net utility
        net_utility = sapply(1:length(genes), function(i) {
            sapply(delta_10_vect, function(d_10) {
               prev[i] * (- d_01 + (d_01 + d_10) * pen[i]) + K
            })
        })
        net_utility = cbind(net_utility, rowSums(net_utility))
        
        # Lower bound of 95% credible interval
        net_utility_lo = sapply(1:length(genes), function(i) {
            sapply(delta_10_vect, function(d_10) {
                p = qbeta(0.025, 
                          shape1 = pen[i] * precision[i], 
                          shape2 = precision[i] - pen[i] * precision[i])
                prev[i] * (- d_01 + (d_01 + d_10) * p) + K
            })
        })
        net_utility_lo = cbind(net_utility_lo, rowSums(net_utility_lo))
        
        # Upper bound of 95% credible interval
        net_utility_hi = sapply(1:length(genes), function(i) {
            sapply(delta_10_vect, function(d_10) {
                p = qbeta(0.975, 
                          shape1 = pen[i] * precision[i], 
                          shape2 = precision[i] - pen[i] * precision[i])
                prev[i] * (- d_01 + (d_01 + d_10) * p) + K
            })
        })
        net_utility_hi = cbind(net_utility_hi, rowSums(net_utility_hi))
        
        # Probability of net utility being positive
        net_prob_pos = sapply(1:length(genes), function(i) {
            sapply(delta_10_vect, function(d_10) {
                1 - pbeta(d_01 / (d_01 + d_10), 
                          shape1 = pen[i] * precision[i], 
                          shape2 = precision[i] - pen[i] * precision[i])
            })
        })
        net_prob_pos = cbind(net_prob_pos, 
                             apply(replicate(100, 
                                             sapply(delta_10_vect, function(d_10) {
                                                 sum(sapply(1:length(genes), function(i) {
                                                     rbeta(1, shape1 = pen[i] * precision[i], 
                                                           shape2 = precision[i] - pen[i] * precision[i])
                                                 }) * prev) / sum(prev) > d_01 / (d_01 + d_10)
                                             })), 1, mean)
        )
        
        # Data for plotting
        thresh_plot_df = 
            do.call(rbind, lapply(1:ncol(net_utility), function(i) {
                data.frame(delta_10 = delta_10_vect, 
                           net_utility = net_utility[,i], 
                           lo = net_utility_lo[,i], 
                           hi = net_utility_hi[,i], 
                           prob_pos = net_prob_pos[,i], 
                           gene = c(genes, "All")[i])
            }))
        
        # Re-scaling factor for secondary axis (probability of positive utility)
        scale_fact = 2 * max(net_utility_hi)
        p = ggplot(thresh_plot_df) + 
            geom_line(aes(x=delta_10, y=net_utility, group = 1,
                          text = paste("FN/FP Utility Cost:", round(delta_10, 5), 
                                       "<br>Net Utility:", round(net_utility, 5)), 
                          color = "Net Utility", linetype = "Net Utility")) + 
            geom_ribbon(aes(x=delta_10, ymin=lo, ymax=hi, group = 1,
                            text = paste("FN/FP Utility Cost:", round(delta_10, 5),  
                                         "<br>2.5% Net Utility:", round(lo, 5), 
                                         "<br>97.5% Net Utility:", round(hi, 5))), 
                        alpha=0.2) + 
            facet_wrap(~gene, ncol = 3) + 
            geom_line(aes(x=delta_10, y = (prob_pos - 0.5) * scale_fact, group = 1, 
                          text = paste("FN/FP Utility Cost:", round(delta_10, 5),  
                                       "<br>Probability of Positive Net Utility:", 
                                       round((prob_pos - 0.5) * scale_fact, 5)), 
                          color = "Probability of Positive Net Utility", 
                          linetype = "Probability of Positive Net Utility")) +
            coord_cartesian(ylim = c(-0.01, 0.05)) + 
            geom_hline(aes(yintercept = 0, 
                           color = "Reference Line", 
                           linetype = "Reference Line"), show.legend = TRUE) + 
            geom_vline(data = data.frame(gene = c(genes, "All"), 
                                         thresh = c(1/pen - 1, 
                                                    sum(prev) / sum(prev*pen) - 1)), 
                       aes(xintercept = thresh, 
                           color = "Reference Line", 
                           linetype = "Reference Line"), show.legend = TRUE) +
            scale_color_manual(values = c("black", "firebrick3", "black"), 
                               labels = c("Net Utility",  
                                          "Probability of Positive Net Utility", 
                                          "Reference Line"), 
                               guide = guide_legend(title = NULL)) +
            scale_linetype_manual(values = c("solid", "dotted", "dashed"), 
                                  labels = c("Net Utility",  
                                             "Probability of Positive Net Utility", 
                                             "Reference Line"), 
                                  guide = guide_legend(title = NULL)) +
            theme(axis.title.y.right = element_text(color = "firebrick3"), 
                  axis.text.y.right = element_text(color = "firebrick3")) +
            xlab("FN/FP Utility Cost") + ylab("Net Utility of Testing for Gene")
        
        # Width of the Shiny app
        width = get_width()
        # Number of rows of plots
        nrows = ceiling((length(genes) + 1)/3)
        p %>% ggplotly(width = width * 0.8, height = width * 0.8 * 2/3, 
                       tooltip = "text") %>% 
            layout(legend = list(orientation = "h", x = 0, y = -0.1))
    })

})

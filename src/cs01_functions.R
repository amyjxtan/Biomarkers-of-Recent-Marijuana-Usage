compound_scatterplot_group <- function(dataset, compound, timepoints){
  if(max(dataset[,compound],na.rm=TRUE)==0){
    print(
      dataset |> 
        filter(!is.na(time_from_start)) |>
        ggplot(aes_string(x="time_from_start", 
                          y=compound,
                          color="group")) + 
        geom_point() +
        geom_vline(data=timepoints, aes(xintercept=as.numeric(stop)), 
                   linetype="dashed", 
                   color="gray28") +
        scale_color_manual(values=c("#19831C", "#A27FC9")) +
        scale_y_continuous(limits=c(0,3)) +
        theme_classic() +
        theme(legend.position="bottom",
              legend.title=element_blank()) +
        labs(x='Time From Start (min)',
             y=gsub('GLUC', 'gluc',gsub("_", "-", toupper(compound))))
    )}else{
      print(
        dataset |> 
          filter(!is.na(time_from_start)) |>
          ggplot(aes_string(x="time_from_start", 
                            y=compound,
                            color="group")) + 
          geom_point() +
          geom_vline(data=timepoints, aes(xintercept=as.numeric(stop)), 
                     linetype="dashed", 
                     color="gray28")  +
          scale_color_manual(values=c("#19831C", "#A27FC9")) +
          theme_classic() +
          theme(legend.position="bottom",
                legend.title=element_blank()) +
          labs(x='Time From Start (min)',
               y=gsub('GLUC', 'gluc', gsub("_", "-", toupper(compound))))
      )
    }
}

make_calculations <- function(dataset, dataset_removedups, split, compound, 
                              start = start, stop = stop, tpt_use = tpt_use){
  ## remove NAs
  df <- dataset_removedups %>% 
    dplyr::select(treatment, compound, timepoint_use) %>%
    rename(compound = 2) %>%
    filter(complete.cases(.))
  if(nrow(df)>0){
    if(stop <= 0){
      output <- df %>% 
        summarise(TP = 0,
                  FN = 0,
                  FP = sum(compound >= split),
                  TN = sum(compound < split)) 
    }else{
      if(split == 0){
        output_pre <- df %>% 
          filter(tpt_use == "pre-smoking") %>%
          summarise(TP = 0,
                    FN = 0,
                    FP = sum(compound >= split),
                    TN = sum(compound < split)) 
        
        output <- df %>% 
          filter(tpt_use != "pre-smoking") %>%
          summarise(TP = sum(treatment != "Placebo" & compound > split),
                    FN = sum(treatment != "Placebo" & compound <= split),
                    FP = sum(treatment == "Placebo" & compound > split),
                    TN = sum(treatment == "Placebo" & compound < split))
        
        output <- output + output_pre
      }else{
        ## calculate values if pre-smoking
        output_pre <- df %>% 
          filter(tpt_use == "pre-smoking") %>%
          summarise(TP = 0,
                    FN = 0,
                    FP = sum(compound >= split),
                    TN = sum(compound < split)) 
        
        output <- df %>% 
          filter(tpt_use != "pre-smoking") %>%
          summarise(TP = sum(treatment != "Placebo" & compound >= split),
                    FN = sum(treatment != "Placebo" & compound < split),
                    FP = sum(treatment == "Placebo" & compound >= split),
                    TN = sum(treatment == "Placebo" & compound < split))
        
        output <- output + output_pre
      }
    }
  }
  # clean things up; make calculations on above values
  output <- output %>%
    mutate(detection_limit = split,
           compound = compound,
           time_start = start,
           time_stop = stop,
           time_window = tpt_use,
           NAs = nrow(dataset) - nrow(df),
           N = nrow(dataset_removedups),
           N_removed = nrow(dataset) - nrow(dataset_removedups),
           Sensitivity = (TP/(TP + FN)), 
           Specificity = (TN /(TN + FP)),
           PPV = (TP/(TP+FP)),
           NPV = (TN/(TN + FN)),
           Efficiency = ((TP + TN)/(TP + TN + FP + FN))*100
    )
  
  return(output)
}

sens_spec <- function(dataset, compound, start, stop, tpt_use, 
                      lowest_value = 0.5, splits = NULL, ...){
  # if it's not all NAs...
  if(sum(is.na(dataset[,compound])) != nrow(dataset)){
    # specify what splits should be used for calculations
    if(is.null(splits)){
      limits <- dataset[is.finite(rowSums(dataset[,compound])),compound]
      ## define lower and upper limits
      lower = min(limits, na.rm=TRUE)
      upper = max(limits, na.rm=TRUE)
      ## determine splits to use for calculations
      tosplit = pull(limits[,1])[limits[,1]>0]
      ## only split if there are detectable limits:
      if(length(tosplit)>=1){
        splits = c(lowest_value, quantile(tosplit, probs=seq(0, 1, by = 0.01), na.rm=TRUE))
        splits = unique(splits)
      }else{
        splits = 0
      }
    }else{
      splits = splits
    }
    # filter to include timepoint of interest
    dataset <- dataset %>% 
      filter(time_from_start > start & time_from_start <= stop & !is.na(timepoint_use))
    dataset_removedups <- dataset %>%
      filter(!is.na(timepoint_use)) %>% 
      group_by(timepoint_use) %>% 
      distinct(id, .keep_all = TRUE) %>% 
      ungroup()

    ## create empty output variable which we'll fill in
    ## iterate through each possible dose and calculate
    output <- map_dfr(as.list(splits), ~make_calculations(dataset, 
                                                          dataset_removedups, 
                                                          split = .x,
                                                          compound,
                                                          start = start,
                                                          stop = stop, 
                                                          tpt_use = tpt_use))
  }
  
  return(output)
}

sens_spec_cpd <- function(dataset, cpd, timepoints, splits = NULL){
  args2 <- list(start = timepoints$start, 
                stop = timepoints$stop, 
                tpt_use = timepoints$timepoint)
  out <- args2 %>% 
    pmap_dfr(sens_spec, dataset, compound = cpd, splits = splits)
  return(out)
}

# helper function to clean up name of two compounds
clean_gluc <- function(df){
  df <- df |> 
    mutate(compound=gsub('GLUC', 'gluc',gsub("_","-",toupper(compound))),
           compound=gsub('THCOH', '11-OH-THC', compound))
  return(df)
}

ss_plot <- function(output, tpts=8, tissue){
  to_include = output %>%
    group_by(compound) %>% 
    summarize(mean_detection = mean(detection_limit)) %>% 
    filter(mean_detection > 0)
  
  output <-  output %>% 
    mutate(iszero = ifelse(time_start<0,TRUE,FALSE),
           Sensitivity = round(Sensitivity*100,0),
           Specificity = round(Specificity*100,0)) %>%
    filter(compound %in% to_include$compound,
           time_window != "pre-smoking") %>%
    clean_gluc() %>% 
    mutate(compound = fct_relevel(as.factor(compound), "THC"))
  
  output <- output %>%  mutate(
    legend = paste0(time_window,' (N=', N,')'))
  
  blue_colors = c('#C2F8FF', '#A2DDED', '#86BEDC', '#6C9FCA', 
                  '#547EB9', '#3F5EA8', '#2D4096', '#1E2385',
                  '#181173', '#180762', '#180051')
  values = c(blue_colors[1:tpts])
  
  print(ggplot(output, aes(x = detection_limit, y = Sensitivity, group = fct_inorder(legend))) + 
          geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
          geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
          facet_grid(~compound, scales = "free_x") +
          labs(x = 'Detection Limit',
               y = 'Sensitivity') +
          ylim(0,1) +
          scale_color_manual(values = values, name = 'Time Window') +
          theme_classic(base_size = 12) + 
          theme(axis.title = element_text(size=16), 
                panel.grid = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_text(size = 12))  
  )
  print(
    ggplot(output, aes(x = detection_limit, y = Specificity, group = fct_inorder(legend))) + 
      geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
      geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
      facet_grid(~compound, scales = "free_x") +
      ylim(0,100) +
      labs(title = tissue,
           x = 'Detection Limit',
           y = 'Specificity') +
      scale_color_manual(values = values, name = 'Time Window') +
      theme_classic(base_size = 12) + 
      theme(axis.title = element_text(size=16),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12))
  )
  print(
    ggplot(output, aes(x=(100-Specificity), y = Sensitivity, group = fct_inorder(legend))) +
      geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
      geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
      facet_grid(~compound) +
      xlim(0, 100) +
      ylim(0, 100) +
      labs(title = tissue,
           x = '(100-Specificity)',
           y = 'Sensitivity') +
      scale_color_manual(values = values, name = 'Time Window') +
      theme_classic(base_size = 12) + 
      theme(axis.title = element_text(size=16),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12),
            axis.text = element_text(size=12))
  )
}

roc_plot <- function(output, tpts=8, tissue){
  to_include = output %>%
    group_by(compound) %>% 
    summarize(mean_detection = mean(detection_limit)) %>% 
    filter(mean_detection > 0)
  
  output <-  output %>% 
    mutate(iszero = ifelse(time_start<0,TRUE,FALSE),
           Sensitivity = round(Sensitivity*100,0),
           Specificity = round(Specificity*100,0)) %>%
    filter(compound %in% to_include$compound,
           time_window != "pre-smoking") %>%
    clean_gluc() %>% 
    mutate(compound = fct_relevel(as.factor(compound), "THC"))
  
  output <- output %>% mutate(
    legend = paste0(time_window,' (N=', N,')'))
  
  blue_colors = c('#C2F8FF', '#86BEDC', 
                  '#547EB9', '#2D4096',
                  '#181173', '#180051')
  values = c(blue_colors[1:tpts])
  print(
    ggplot(output, aes(x=(100-Specificity), y = Sensitivity, group = fct_inorder(legend))) +
      geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
      geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
      facet_grid(~compound) +
      xlim(0, 100) +
      ylim(0, 100) +
      labs(title = tissue,
           x = '(100-Specificity)',
           y = 'Sensitivity') +
      scale_color_manual(values = values, name = 'Time Window') +
      theme_classic(base_size = 12) + 
      theme(axis.title = element_text(size=16),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12),
            axis.text = element_text(size=12) )
  )
}
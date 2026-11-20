# format(): the text style is what print() shows; digits; unknown style

    Code
      print(et, digits = 2)
    Output
       response type        parameter                mean (sd)      95% interval  
       y        fixed       (Intercept)              -1.1 (0.1)     [-1.3, -0.87] 
                fixed       t                        -0.023 (0.073) [-0.17, 0.12] 
                association gamma (x1 -> y)          0.91 (0.081)   [0.75, 1.1]   
                family      phi                      23 (2.7)       [19, 29]      
       x1       fixed       (Intercept)              0.99 (0.081)   [0.83, 1.2]   
                fixed       t                        0.47 (0.049)   [0.37, 0.57]  
                sd          SD (Intercept) | id      0.49 (0.057)   [0.38, 0.61]  
                sd          SD t | id                0.26 (0.034)   [0.2, 0.34]   
                correlation corr (Intercept), t | id -0.27 (0.16)   [-0.57, 0.062]
                residual    residual SD              0.18 (0.013)   [0.15, 0.21]  

# format(style = 'markdown'): a pipe table, `|` of labels escaped, response on its first row

    Code
      cat(format(et, "markdown", digits = 2))
    Output
      | response | parameter | mean (sd) | 95% interval |
      |---|---|---|---|
      | y | (Intercept) | -1.1 (0.1) | [-1.3, -0.87] |
      |  | t | -0.023 (0.073) | [-0.17, 0.12] |
      |  | gamma (x1 -> y) | 0.91 (0.081) | [0.75, 1.1] |
      |  | phi | 23 (2.7) | [19, 29] |
      | x1 | (Intercept) | 0.99 (0.081) | [0.83, 1.2] |
      |  | t | 0.47 (0.049) | [0.37, 0.57] |
      |  | SD (Intercept) \| id | 0.49 (0.057) | [0.38, 0.61] |
      |  | SD t \| id | 0.26 (0.034) | [0.2, 0.34] |
      |  | corr (Intercept), t \| id | -0.27 (0.16) | [-0.57, 0.062] |
      |  | residual SD | 0.18 (0.013) | [0.15, 0.21] |

# format(style = 'latex'): a tabular, -> as $\to$, special characters escaped

    Code
      cat(format(et, "latex", digits = 2))
    Output
      \begin{tabular}{llrr}
      \hline
      response & parameter & mean (sd) & 95\% interval \\
      \hline
      y & (Intercept) & -1.1 (0.1) & [-1.3, -0.87] \\
       & t & -0.023 (0.073) & [-0.17, 0.12] \\
       & gamma (x1 $\to$ y) & 0.91 (0.081) & [0.75, 1.1] \\
       & phi & 23 (2.7) & [19, 29] \\
      x1 & (Intercept) & 0.99 (0.081) & [0.83, 1.2] \\
       & t & 0.47 (0.049) & [0.37, 0.57] \\
       & SD (Intercept) | id & 0.49 (0.057) & [0.38, 0.61] \\
       & SD t | id & 0.26 (0.034) & [0.2, 0.34] \\
       & corr (Intercept), t | id & -0.27 (0.16) & [-0.57, 0.062] \\
       & residual SD & 0.18 (0.013) & [0.15, 0.21] \\
      \hline
      \end{tabular}

# format(style = 'array'): the table as display math, cells in \text{}

    Code
      cat(format(et, "array", digits = 2))
    Output
      \begin{array}{llrr}
      \hline
      \text{response} & \text{parameter} & \text{mean (sd)} & \text{95\% interval} \\
      \hline
      \text{y} & \text{(Intercept)} & -1.1\ (0.1) & [-1.3,\ -0.87] \\
       & \text{t} & -0.023\ (0.073) & [-0.17,\ 0.12] \\
       & \text{gamma (x1 $\to$ y)} & 0.91\ (0.081) & [0.75,\ 1.1] \\
       & \text{phi} & 23\ (2.7) & [19,\ 29] \\
      \text{x1} & \text{(Intercept)} & 0.99\ (0.081) & [0.83,\ 1.2] \\
       & \text{t} & 0.47\ (0.049) & [0.37,\ 0.57] \\
       & \text{SD (Intercept) | id} & 0.49\ (0.057) & [0.38,\ 0.61] \\
       & \text{SD t | id} & 0.26\ (0.034) & [0.2,\ 0.34] \\
       & \text{corr (Intercept), t | id} & -0.27\ (0.16) & [-0.57,\ 0.062] \\
       & \text{residual SD} & 0.18\ (0.013) & [0.15,\ 0.21] \\
      \hline
      \end{array}


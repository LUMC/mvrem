# format(): the text style is what print() shows; digits; unknown style

    Code
      cat(maskNumbers(format(et, "text", digits = 2)))
    Output
       response type parameter mean (sd) #% interval 
       y fixed (Intercept) # (#) [#, #] 
       fixed t # (#) [#, #] 
       association gamma (x1 -> y) # (#) [#, #] 
       family phi # (#) [#, #] 
       x1 fixed (Intercept) # (#) [#, #] 
       fixed t # (#) [#, #] 
       sd SD (Intercept) | id # (#) [#, #] 
       sd SD t | id # (#) [#, #] 
       correlation corr (Intercept), t | id # (#) [#, #]
       residual residual SD # (#) [#, #] 

# format(style = 'markdown'): a pipe table, `|` of labels escaped, response on its first row

    Code
      cat(maskNumbers(format(et, "markdown", digits = 2)))
    Output
      | response | parameter | mean (sd) | #% interval |
      |---|---|---|---|
      | y | (Intercept) | # (#) | [#, #] |
      | | t | # (#) | [#, #] |
      | | gamma (x1 -> y) | # (#) | [#, #] |
      | | phi | # (#) | [#, #] |
      | x1 | (Intercept) | # (#) | [#, #] |
      | | t | # (#) | [#, #] |
      | | SD (Intercept) \| id | # (#) | [#, #] |
      | | SD t \| id | # (#) | [#, #] |
      | | corr (Intercept), t \| id | # (#) | [#, #] |
      | | residual SD | # (#) | [#, #] |

# format(style = 'latex'): a tabular, -> as $\to$, special characters escaped

    Code
      cat(maskNumbers(format(et, "latex", digits = 2)))
    Output
      \begin{tabular}{llrr}
      \hline
      response & parameter & mean (sd) & #\% interval \\
      \hline
      y & (Intercept) & # (#) & [#, #] \\
       & t & # (#) & [#, #] \\
       & gamma (x1 $\to$ y) & # (#) & [#, #] \\
       & phi & # (#) & [#, #] \\
      x1 & (Intercept) & # (#) & [#, #] \\
       & t & # (#) & [#, #] \\
       & SD (Intercept) | id & # (#) & [#, #] \\
       & SD t | id & # (#) & [#, #] \\
       & corr (Intercept), t | id & # (#) & [#, #] \\
       & residual SD & # (#) & [#, #] \\
      \hline
      \end{tabular}

# format(style = 'array'): the table as display math, cells in \text{}

    Code
      cat(maskNumbers(format(et, "array", digits = 2)))
    Output
      \begin{array}{llrr}
      \hline
      \text{response} & \text{parameter} & \text{mean (sd)} & \text{#\% interval} \\
      \hline
      \text{y} & \text{(Intercept)} & #\ (#) & [#,\ #] \\
       & \text{t} & #\ (#) & [#,\ #] \\
       & \text{gamma (x1 $\to$ y)} & #\ (#) & [#,\ #] \\
       & \text{phi} & #\ (#) & [#,\ #] \\
      \text{x1} & \text{(Intercept)} & #\ (#) & [#,\ #] \\
       & \text{t} & #\ (#) & [#,\ #] \\
       & \text{SD (Intercept) | id} & #\ (#) & [#,\ #] \\
       & \text{SD t | id} & #\ (#) & [#,\ #] \\
       & \text{corr (Intercept), t | id} & #\ (#) & [#,\ #] \\
       & \text{residual SD} & #\ (#) & [#,\ #] \\
      \hline
      \end{array}


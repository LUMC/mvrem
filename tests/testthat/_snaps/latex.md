# latexEquation(): snapshot of the whole equation of each fixture

    Code
      cat(format(latexEquation(fxModel(nm, fit = FALSE))))
    Output
      \begin{aligned}
      x_{1,i} &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t_{i} + b^{(x_{1})}_{0,\mathrm{id}[i]} + b^{(x_{1})}_{1,\mathrm{id}[i]}\, t_{i} + \varepsilon^{(x_{1})}_{i} \\
      0 &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t + b^{(x_{1})}_{0,\mathrm{id}} + b^{(x_{1})}_{1,\mathrm{id}}\, t - u^{(x_{1} \to y)}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{y} \\
      \operatorname{logit}\left(\mu^{(y)}_{i}\right) &= \beta^{(y)}_{0} + \beta^{(y)}_{1}\, t_{i} + \gamma_{y,x_{1}}\,u^{(x_{1} \to y)}_{\mathrm{id}[i],t_{i}} \\
      \varepsilon^{(x_{1})}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{1})}_{\varepsilon}\right)^{-1}\right) \\
      \left(b^{(x_{1})}_{0,l}, b^{(x_{1})}_{1,l}\right)^\top &\sim \mathcal{N}\left(\mathbf{0}, \left(\mathbf{Q}^{(x_{1})}_{\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      y_{i} &\sim \mathrm{Beta}\left(\mu^{(y)}_{i}\, \phi^{(y)}, \left(1 - \mu^{(y)}_{i}\right) \phi^{(y)}\right) \\
      u^{(x_{1} \to y)}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{1} \to y}\right)^{-1}\right)
      \end{aligned}

---

    Code
      cat(format(latexEquation(fxModel(nm, fit = FALSE))))
    Output
      \begin{aligned}
      x_{1,i} &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t_{i} + b^{(x_{1})}_{0,\mathrm{id}[i]} + b^{(x_{1})}_{1,\mathrm{id}[i]}\, t_{i} + \varepsilon^{(x_{1})}_{i} \\
      \log\left(\lambda^{(x_{2})}_{i}\right) &= \beta^{(x_{2})}_{0} + \beta^{(x_{2})}_{1}\, t_{i} + b^{(x_{2})}_{0,\mathrm{id}[i]} \\
      0 &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t + b^{(x_{1})}_{0,\mathrm{id}} + b^{(x_{1})}_{1,\mathrm{id}}\, t - u^{(x_{1} \to k)}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{k} \\
      0 &= \beta^{(x_{2})}_{0} + \beta^{(x_{2})}_{1}\, t + b^{(x_{2})}_{0,\mathrm{id}} - u^{(x_{2} \to k)}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{k} \\
      0 &= \beta^{(x_{2})}_{0} + \beta^{(x_{2})}_{1}\, t + b^{(x_{2})}_{0,\mathrm{id}} - u^{(x_{2} \to \mathit{cnt})}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{\mathit{cnt}} \\
      \operatorname{logit}\left(p^{(k)}_{i}\right) &= \beta^{(k)}_{0} + \beta^{(k)}_{1}\, t_{i} + \gamma_{k,x_{1}}\,u^{(x_{1} \to k)}_{\mathrm{id}[i],t_{i}} + \gamma_{k,x_{2}}\,u^{(x_{2} \to k)}_{\mathrm{id}[i],t_{i}} \\
      \log\left(\lambda^{(\mathit{cnt})}_{i}\right) &= \beta^{(\mathit{cnt})}_{0} + \beta^{(\mathit{cnt})}_{1}\, \mathbb{I}\left[\mathit{grp}_{i} = \text{b}\right] + \beta^{(\mathit{cnt})}_{2}\, \mathbb{I}\left[\mathit{grp}_{i} = \text{c}\right] + \gamma_{\mathit{cnt},x_{2}}\,u^{(x_{2} \to \mathit{cnt})}_{\mathrm{id}[i],t_{i}} \\
      \varepsilon^{(x_{1})}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{1})}_{\varepsilon}\right)^{-1}\right) \\
      \left(b^{(x_{1})}_{0,l}, b^{(x_{1})}_{1,l}\right)^\top &\sim \mathcal{N}\left(\mathbf{0}, \left(\mathbf{Q}^{(x_{1})}_{\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      x_{2,i} &\sim \mathrm{Poisson}\left(\lambda^{(x_{2})}_{i}\right) \\
      b^{(x_{2})}_{0,l} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{2})}_{0,\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      k_{i} &\sim \mathrm{Binomial}\left(n_{i}, p^{(k)}_{i}\right) \\
      \mathit{cnt}_{i} &\sim \mathrm{Poisson}\left(\mathit{pt}_{i}\, \lambda^{(\mathit{cnt})}_{i}\right) \\
      u^{(x_{1} \to k)}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{1} \to k}\right)^{-1}\right) \\
      u^{(x_{2} \to k)}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{2} \to k}\right)^{-1}\right) \\
      u^{(x_{2} \to \mathit{cnt})}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{2} \to \mathit{cnt}}\right)^{-1}\right)
      \end{aligned}

---

    Code
      cat(format(latexEquation(fxModel(nm, fit = FALSE))))
    Output
      \begin{aligned}
      x_{1,i} &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t_{i} + b^{(x_{1})}_{0,\mathrm{id}[i]} + b^{(x_{1})}_{1,\mathrm{id}[i]}\, t_{i} + \varepsilon^{(x_{1})}_{i} \\
      0 &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t + b^{(x_{1})}_{0,\mathrm{id}} + b^{(x_{1})}_{1,\mathrm{id}}\, t - u^{(x_{1} \to \mathit{time})}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{\mathit{time}} \\
      \log\left(\lambda^{(\mathit{time})}_{i}\right) &= \beta^{(\mathit{time})}_{0} + \gamma_{\mathit{time},x_{1}}\,u^{(x_{1} \to \mathit{time})}_{\mathrm{id}[i],t_{i}} \\
      \varepsilon^{(x_{1})}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{1})}_{\varepsilon}\right)^{-1}\right) \\
      \left(b^{(x_{1})}_{0,l}, b^{(x_{1})}_{1,l}\right)^\top &\sim \mathcal{N}\left(\mathbf{0}, \left(\mathbf{Q}^{(x_{1})}_{\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      \mathit{time}_{i} &\sim \mathrm{Weibull}\left(\alpha^{(\mathit{time})}, \lambda^{(\mathit{time})}_{i}\right), \quad S(t) = \exp\left(-\lambda^{(\mathit{time})}_{i}\, t^{\alpha^{(\mathit{time})}}\right), \quad \mathit{status}_{i} \in \{0, 1\} \\
      u^{(x_{1} \to \mathit{time})}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{1} \to \mathit{time}}\right)^{-1}\right)
      \end{aligned}

---

    Code
      cat(format(latexEquation(fxModel(nm, fit = FALSE))))
    Output
      \begin{aligned}
      x_{1,i} &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t_{i} + b^{(x_{1})}_{0,\mathrm{id}[i]} + b^{(x_{1})}_{1,\mathrm{id}[i]}\, t_{i} + \varepsilon^{(x_{1})}_{i} \\
      \log\left(\mu^{(x_{2})}_{i}\right) &= \beta^{(x_{2})}_{0} + \beta^{(x_{2})}_{1}\, t_{i} + b^{(x_{2})}_{0,\mathrm{id}[i]} + b^{(x_{2})}_{1,\mathrm{id}[i]}\, t_{i} \\
      x_{3,i} &= \beta^{(x_{3})}_{0} + \beta^{(x_{3})}_{1}\, t_{i} + b^{(x_{3})}_{0,\mathrm{id}[i]} + b^{(x_{3})}_{1,\mathrm{id}[i]}\, t_{i} + \varepsilon^{(x_{3})}_{i} \\
      0 &= \beta^{(x_{1})}_{0} + \beta^{(x_{1})}_{1}\, t + b^{(x_{1})}_{0,\mathrm{id}} + b^{(x_{1})}_{1,\mathrm{id}}\, t - u^{(x_{1} \to \mathit{yg})}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{\mathit{yg}} \\
      0 &= \beta^{(x_{2})}_{0} + \beta^{(x_{2})}_{1}\, t + b^{(x_{2})}_{0,\mathrm{id}} + b^{(x_{2})}_{1,\mathrm{id}}\, t - u^{(x_{2} \to \mathit{yg})}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{\mathit{yg}} \\
      0 &= \beta^{(x_{3})}_{0} + \beta^{(x_{3})}_{1}\, t + b^{(x_{3})}_{0,\mathrm{id}} + b^{(x_{3})}_{1,\mathrm{id}}\, t - u^{(x_{3} \to \mathit{yg})}_{\mathrm{id},t}, \quad \left(\mathrm{id}, t\right) \in \mathcal{D}_{\mathit{yg}} \\
      \mathit{yg}_{i} &= \beta^{(\mathit{yg})}_{0} + \beta^{(\mathit{yg})}_{1}\, t_{i} + b^{(\mathit{yg})}_{0,\mathrm{site}[i]} + \gamma_{\mathit{yg},x_{1}}\,u^{(x_{1} \to \mathit{yg})}_{\mathrm{id}[i],t_{i}} + \gamma_{\mathit{yg},x_{2}}\,u^{(x_{2} \to \mathit{yg})}_{\mathrm{id}[i],t_{i}} + \gamma_{\mathit{yg},x_{3}}\,u^{(x_{3} \to \mathit{yg})}_{\mathrm{id}[i],t_{i}} + \varepsilon^{(\mathit{yg})}_{i} \\
      \varepsilon^{(x_{1})}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{1})}_{\varepsilon}\right)^{-1}\right) \\
      \left(b^{(x_{1})}_{0,l}, b^{(x_{1})}_{1,l}\right)^\top &\sim \mathcal{N}\left(\mathbf{0}, \left(\mathbf{Q}^{(x_{1})}_{\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      \mathbf{Q}^{(x_{1})}_{\mathrm{id}} &\sim \mathcal{W}_{2}\left(4, \begin{pmatrix}1.3 & 1.27 \\ 1.27 & 7.19\end{pmatrix}\right) \\
      x_{2,i} &\sim \mathrm{NegBin}\left(\mu^{(x_{2})}_{i}, \kappa^{(x_{2})}\right) \\
      \kappa^{(x_{2})} &\sim \mathrm{Gamma}\left(1, 0.01\right) \\
      b^{(x_{2})}_{0,l} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{2})}_{0,\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      \left(\tau^{(x_{2})}_{0,\mathrm{id}}\right)^{-1/2} &\sim \mathrm{C}^{+}\left(0, 1\right) \\
      b^{(x_{2})}_{1,l} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{2})}_{1,\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      \left(\tau^{(x_{2})}_{1,\mathrm{id}}\right)^{-1/2} &\sim \mathcal{N}^{+}\left(0, \left(4\right)^{-1}\right) \\
      \varepsilon^{(x_{3})}_{i} &\sim t_{\nu^{(x_{3})}}\left(0, \left(\tau^{(x_{3})}_{\varepsilon}\right)^{-1}\right) \\
      \tau^{(x_{3})}_{\varepsilon} &\sim \mathrm{Gamma}\left(1, 5 \cdot 10^{-5}\right) \\
      b^{(x_{3})}_{0,l} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{3})}_{0,\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      b^{(x_{3})}_{1,l} &\sim \mathcal{N}\left(0, \left(\tau^{(x_{3})}_{1,\mathrm{id}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{id}} \\
      \tau^{(x_{3})}_{k,\mathrm{id}} &\sim \mathrm{Gamma}\left(1, 0.01\right), \quad k = 0, 1 \\
      \varepsilon^{(\mathit{yg})}_{i} &\sim \mathcal{N}\left(0, \left(\tau^{(\mathit{yg})}_{\varepsilon}\right)^{-1}\right) \\
      \left(\tau^{(\mathit{yg})}_{\varepsilon}\right)^{-1/2} &\sim \mathrm{Exp}\left(4.61\right), \quad P\left(\left(\tau^{(\mathit{yg})}_{\varepsilon}\right)^{-1/2} > 1\right) = 0.01 \\
      b^{(\mathit{yg})}_{0,l} &\sim \mathcal{N}\left(0, \left(\tau^{(\mathit{yg})}_{0,\mathrm{site}}\right)^{-1}\right), \quad l = 1, \ldots, n_{\mathrm{site}} \\
      \left(\tau^{(\mathit{yg})}_{0,\mathrm{site}}\right)^{-1/2} &\sim \mathrm{Exp}\left(4.61\right), \quad P\left(\left(\tau^{(\mathit{yg})}_{0,\mathrm{site}}\right)^{-1/2} > 1\right) = 0.01 \\
      \gamma_{\mathit{yg},x_{1}} &\sim \mathcal{N}\left(0, \left(0.25\right)^{-1}\right) \\
      u^{(x_{1} \to \mathit{yg})}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{1} \to \mathit{yg}}\right)^{-1}\right) \\
      u^{(x_{2} \to \mathit{yg})}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{2} \to \mathit{yg}}\right)^{-1}\right) \\
      u^{(x_{3} \to \mathit{yg})}_{\mathrm{id},t} &\sim \mathcal{N}\left(0, \left(\tau_{x_{3} \to \mathit{yg}}\right)^{-1}\right)
      \end{aligned}


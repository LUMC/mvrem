# syntax=docker/dockerfile:1
FROM rocker/geospatial:4.5.2

RUN mkdir -p /home/rstudio/project && \ 
    chown -R rstudio:rstudio /home/rstudio/project
WORKDIR /home/rstudio/project

RUN R -e "options(repos = c(INLA = 'https://inla.r-inla-download.org/R/testing', CRAN = 'https://cloud.r-project.org')); install.packages('INLA')" && \
    chmod -R a+x /usr/local/lib/R/site-library/INLA/bin/linux && \
    R -e 'n = 100; a = 1; b = 1; tau = 100; z = rnorm(n); eta = a + b*z; scale = exp(rnorm(n)); prec = scale*tau; y = rnorm(n, mean = eta, sd = 1/sqrt(prec)); result = INLA::inla(y ~ 1+z, family = "gaussian", data = list(y=y, z=z)); summary(result)'

RUN R -e "install.packages(c('languageserver', 'httpgd', 'httpuv', 'jsonlite', 'Rlab', 'lme4', 'reactable'))" && \
    R -e "install.packages(c('jgd'), repos = c('https://reditorsupport.r-universe.dev', getOption('repos')))"

# Install arf CLI
RUN curl --proto '=https' --tlsv1.2 -LsSf https://github.com/eitsupi/arf/releases/download/v0.5.2/arf-console-installer.sh | sh

WORKDIR /home/rstudio

##COPY <<EOT /home/rstudio/.Rprofile
#if (interactive() && Sys.getenv("RSTUDIO") == "") {
#    Sys.setenv(TERM_PROGRAM = "vscode")
#    source(file.path(Sys.getenv("HOME"), ".vscode-R", "init.R"))
#}
#EOT
#RUN chown rstudio:rstudio /home/rstudio/.Rprofile


#RUN chmod -R u+w /home/rstudio/project && rm -rf /home/rstudio/project

EXPOSE 8787
CMD ["/init"]

# To prepare the Docker image, run the following commands in the top-level repository directory:
#   docker build -t="inla" --platform linux/amd64 containers
#   docker images
#
# To run R console (with working library INLA):
#   docker run --name inla --rm -ti -v .:/home/rstudio inla R
#
# To run RStudio server at localhost:8787:
#   docker run --name inla --rm -ti -v .:/home/rstudio -p 8787:8787 inla




# To prepare the dockercontainer image, run the following commands in the top-level repository directory:
#   #cd GIT-MAIN-DIRECTORY    
#   devcontainer build --platform=linux/amd64 --workspace-folder .  >.devcontainer/dc.json    # --no-cache
#   IMAGE_ID=$(jq -r '.imageName[0]' .devcontainer/dc.json)
#   CONTAINER_ID=$(docker run --rm -d -v .:/home/rstudio/project ${IMAGE_ID})
#   code --folder-uri vscode-remote://attached-container+$(printf "$CONTAINER_ID" | od -A n -t x1 | sed 's/ *//g' | tr -d '\n')/home/rstudio/project
#   docker container kill ${CONTAINER_ID}

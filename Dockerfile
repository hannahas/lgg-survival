# rocker/shiny is a purpose-built, actively maintained image for
# R + Shiny apps, built on a stable Debian release rather than
# Debian's rolling testing/sid repos — avoids the package-conflict
# issue that r-base:4.3.3 hit above. It also ships with Shiny
# Server already configured.
FROM --platform=linux/amd64 rocker/shiny:4.3

# System-level libraries needed to compile survminer/ggplot2's
# dependencies. Fewer are needed here than before, since
# rocker/shiny already includes several common ones.
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    gfortran \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install the R packages the app needs
RUN R -e "install.packages(c('shiny', 'survival', 'ggplot2'), repos='https://cloud.r-project.org'); \
    if (!all(c('shiny','survival','ggplot2') %in% rownames(installed.packages()))) quit(status=1)"

# rocker/shiny expects apps in this specific directory, and its
# built-in Shiny Server automatically serves whatever is placed
# here — replacing the app's own runApp() call below.
COPY app/ /srv/shiny-server/app/

EXPOSE 3838

# rocker/shiny's default CMD already starts Shiny Server, which
# auto-serves everything under /srv/shiny-server/. No custom CMD
# needed — but explicitly running the app also works if you'd
# rather keep it simple and match what you tested locally:
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/app', host='0.0.0.0', port=3838)"]
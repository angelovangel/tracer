# Use a base image with broader platform support, pinning to a version is good practice
FROM --platform=linux/amd64 rocker/shiny-verse:4.5.0

# Copy the lockfile to restore packages
COPY renv.lock .

# Restore the R environment using the lockfile
# This ensures the same package versions are installed
RUN R -e "install.packages('renv')"
RUN R -e "renv::restore()"

# Copy the Shiny app file to the server directory
COPY app.R /srv/shiny-server/
COPY global.R /srv/shiny-server/


# Expose the application port
EXPOSE 3838
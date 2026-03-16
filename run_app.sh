docker run -it \
    -p 3838:3838 \
    --name pardal-appx \
    -v $PWD:/srv/shiny-server/pardal \
    -e SHINY_LOG_LEVEL=TRACE \
    pardal-app


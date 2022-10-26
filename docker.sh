docker pull r-base
docker build -t cmmr/rbiom-base --no-cache - >build_log.txt <<-"EOF"
	FROM r-base
	MAINTAINER Daniel Smith "dansmith@orst.edu"
 
	RUN                                                                           \
	                                                                              \
	  apt-get update                                                              \
	                                                                              \
	  && apt-get install -y                                                       \
	      libudunits2-dev libssl-dev libxml2-dev libcurl4-openssl-dev libgdal-dev \
	                                                                              \
	  && su - -c "R -e \"install.packages(dependencies=TRUE, c(                   \
	      'ape', 'broom', 'ggbeeswarm', 'ggdensity', 'ggnewscale', 'ggpattern',   \
	      'ggplot2', 'ggrepel', 'ggtext', 'glue', 'gridpattern', 'jsonlite',      \
	      'labeling', 'magrittr', 'memoise', 'methods', 'openxlsx', 'optparse',   \
	      'paletteer', 'patchwork', 'plyr', 'Rcpp', 'RcppParallel',               \ 
	      'rlang', 'R.utils', 'slam', 'tidyr', 'tsne', 'uwot', 'vegan',           \
	      'BiocManager', 'remotes' ))\""                                          \
	                                                                              \
	  && su - -c "R -e \"BiocManager::install(c('ggtree', 'rhdf5'))\""            \
	                                                                              \
	  && rm -rf /tmp/* /var/lib/apt/lists/*
EOF



docker build -t cmmr/rbiom --no-cache - >build_log.txt 2>&1 <<-"EOF"
	FROM cmmr/rbiom-base
	MAINTAINER Daniel Smith "dansmith@orst.edu"
 
 
	 RUN su - -c "R -e \"                                                  \
	  remotes::install_github('cmmr/rbiom', upgrade = 'never') \""         \
	                                                                       \
	  && rm -rf /tmp/*
EOF



docker run --rm -it cmmr/rbiom /bin/bash


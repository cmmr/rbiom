docker run -it --rm -v /gpfs1:/gpfs1 cmmr/rbiom /bin/bash
mkdir /cmmr
ln -s /usr/bin /cmmr/bin
ln -s /usr/local/lib/R/site-library/wrapt/exec/deliver.r /usr/bin
/usr/local/lib/R/site-library/rbiom




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
	      pandoc librsvg2-dev libcairo2-dev libharfbuzz-dev libfribidi-dev        \
	      libxt-dev	                                                              \
	                                                                              \
	  && R -e "install.packages(dependencies=TRUE, c(                             \
	      'ape', 'broom', 'ggbeeswarm', 'ggdensity', 'ggnewscale', 'ggpattern',   \
	      'ggplot2', 'ggrepel', 'ggtext', 'glue', 'gridpattern', 'jsonlite',      \
	      'labeling', 'magrittr', 'memoise', 'methods', 'openxlsx', 'optparse',   \
	      'paletteer', 'patchwork', 'plyr', 'Rcpp', 'RcppParallel',               \
	      'rlang', 'R.utils', 'slam', 'tidyr', 'tsne', 'uwot', 'vegan',           \
	      'BiocManager', 'remotes', 'flexdashboard', 'plotly', 'ggdendro' ))"     \
	                                                                              \
	  && R -e "BiocManager::install(c('ggtree', 'rhdf5', 'treeio'))"              \
	                                                                              \
	  && R -e "remotes::install_github('AnalytixWare/ShinySky')"                  \
	                                                                              \
	  && wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip            \
	  && unzip awscli-exe-linux-x86_64.zip                                        \
	  && ./aws/install                                                            \
	                                                                              \
	  && rm -rf aws* /tmp/* /var/lib/apt/lists/*
EOF



docker build -t cmmr/rbiom --no-cache - >build_log.txt 2>&1 <<-"EOF"
	FROM cmmr/rbiom-base
	MAINTAINER Daniel Smith "dansmith@orst.edu"
 
 
	 RUN R -e "remotes::install_github('cmmr/rbiom', upgrade = 'never')"    \
	                                                                        \
	  && R -e "remotes::install_github('cmmr/wrapt', upgrade = 'never',     \
	    auth_token = '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')"           \
	                                                                        \
	  && mkdir -p /cmmr/bin                                                 \
	  && ln -s /usr/bin/Rscript /cmmr/bin                                   \
	  && ln -s /usr/local/lib/R/site-library/wrapt/exec/deliver.r /usr/bin  \
	                                                                        \
	  && rm -rf /tmp/*
EOF


docker run --rm -it cmmr/rbiom /bin/bash


scp -rp /gpfs1/projects/Pools/BiT_requests/MeyerK-PQ00381-20221212 dpsmith@cmmr-web2:/gpfs1/projects/Pools/BiT_requests
scp -rp /gpfs1/projects/Pools/16SV4/Pool1568/Kommagani_Ramakrishna_PQ00295/amplicon_output dpsmith@cmmr-web2:/gpfs1/projects/Pools/16SV4/Pool1568/Kommagani_Ramakrishna_PQ00295

docker run -it --rm -v /gpfs1:/gpfs1 cmmr/rbiom /bin/bash


deliver.r -d /gpfs1/projects/Pools/BiT_requests/RedondoM-PQ00230-20221216/amplicon_output/uparse          -t "PQ00230" -i "Maria Redondo"            -m "Juwan Cormier" -e juwan.cormier@bcm.edu -a uparse -p V4
deliver.r -d /gpfs1/projects/Pools/BiT_requests/VonMoltkeJ-PQ00256-20221214/amplicon_output/uparse        -t "PQ00256" -i "Jakob von Moltke"         -m "Anaid Reyes"   -e agreyes@bcm.edu       -a uparse -p V4
deliver.r -d /gpfs1/projects/Pools/16SV4/Pool1565/Anandasabapathy_Sharmila_PQ00379/amplicon_output/uparse -t "PQ00379" -i "Sharmila Anandasabapathy" -m "Matt Ross"     -e mcross@bcm.edu        -a uparse -p V4
deliver.r -d /gpfs1/projects/Pools/16SV4/Pool1565/Klopp_Ann_PQ00409/amplicon_output/uparse                -t "PQ00409" -i "Ann Klopp"                -m "Anaid Reyes"   -e agreyes@bcm.edu       -a uparse -p V4
deliver.r -d /gpfs1/projects/Pools/16SV4/Pool1565/Patras_Katy_PQ00398/amplicon_output/uparse              -t "PQ00398" -i "Katy Patras"              -m "Anaid Reyes"   -e agreyes@bcm.edu       -a uparse -p V4
deliver.r -d /gpfs1/projects/Pools/16SV4/Pool1565/Riquelme_Erick_PQ00302/amplicon_output/uparse           -t "PQ00302" -i "Erick Riquelme"           -m "Anaid Reyes"   -e agreyes@bcm.edu       -a uparse -p V4
deliver.r -d /gpfs1/projects/Pools/16SV4/Pool1573/Dumas_Orianne_PQ00359/amplicon_output/uparse            -t "PQ00359" -i "Orianne Dumas"            -m "Anaid Reyes"   -e agreyes@bcm.edu       -a uparse -p V4
deliver.r -d /gpfs1/projects/Pools/16SV1V3/Pool1571/Villapol_Sonia_PQ00371/amplicon_output/deblur         -t "PQ00371" -i "Sonia Villapol"           -m "Juwan Cormier" -e juwan.cormier@bcm.edu -a deblur -p V4

`PQ00230  Maria Redondo             https://data.jplab.net/2fXQvZID/index.html`
`PQ00256  Jakob von Moltke          https://data.jplab.net/OPVEQ9V5/index.html`
`PQ00379  Sharmila Anandasabapathy  https://data.jplab.net/unqdQi9D/index.html`
`PQ00409  Ann Klopp                 https://data.jplab.net/5nHAcQxU/index.html`
`PQ00398  Katy Patras               https://data.jplab.net/b9yNnQux/index.html`
`PQ00302  Erick Riquelme            https://data.jplab.net/46fio3cz/index.html`
`PQ00359  Orianne Dumas             https://data.jplab.net/dLneGhxZ/index.html`
PQ00371  Sonia Villapol            https://data.jplab.net//index.html


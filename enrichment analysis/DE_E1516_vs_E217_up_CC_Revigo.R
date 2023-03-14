# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
library(ggrepel)
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0030286","dynein complex",0.131930369727061,-1.81008031573268,-2.29480403900491,4.42135709851314,-14.7980993185234,0.37136903352721,0),
c("GO:0097542","ciliary tip",0.00230518118724132,3.48479372993267,2.6401791541045,2.66464197555613,-3.26000683313858,0.508802043498527,0.1157287),
c("GO:0016935","glycine-gated chloride channel complex",0.000550043233398145,-1.11073716850614,-5.76125167871985,2.04532297878666,-3.17830801158322,0.77977868703549,0.19761107),
c("GO:0031298","replication fork protection complex",0.0078506170585008,-3.83241308875572,1.03262932591884,3.19617618503997,-4.04591658460789,0.508402751348257,0.27816506),
c("GO:0005782","peroxisomal matrix",0.010580831653368,-4.39918130444373,4.0679269236901,3.32572085801941,-3.51918509906018,0.675631553516285,0.36832144),
c("GO:0000783","nuclear telomere cap complex",0.0116209134037935,-3.68003036344132,-1.29413330884886,3.36642295722597,-3.05407491554185,0.45699963246638,0.43282425),
c("GO:0032589","neuron projection membrane",0.00973576523114717,5.3670156070804,0.680613507148416,3.2895889525426,-3.13725748555333,0.615046319242864,0.48232973),
c("GO:0036064","ciliary basal body",0.0448035215567944,1.05858109328153,1.02877035514936,3.95235647732379,-5.17208158686321,0.30870734419332,0.62550126),
c("GO:0031907","microbody lumen",0.010580831653368,-3.6894231234742,4.46603549674701,3.32572085801941,-3.51918509906018,0.675631553516285,0.63431391),
c("GO:0000811","GINS complex",0.0153062030675611,-3.36769093182687,0.700886746073015,3.48600518636224,-3.05950477515596,0.459270642802393,0.6390631));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = log_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 1, ];
p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("/path_to_your_file/revigo-plot.pdf");
ggsave("figures/DE_E1516_vs_E217_up_CC_REVIGO.png", width=8, height=8, dpi=1200);

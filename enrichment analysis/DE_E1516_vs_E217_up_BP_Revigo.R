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
revigo.data <- rbind(c("GO:0007423","sensory organ development",0.150496402810338,1.71068759697397,-6.79214564921435,4.5602295339144,-3.70601704334959,0.955684293389268,0),
c("GO:0030203","glycosaminoglycan metabolic process",1.00909565689903,-2.83955529539667,-2.10691784416848,5.38662557915907,-3.6219276590108,0.986009109205803,0),
c("GO:0060271","cilium assembly",0.155248353865325,-5.65977503341237,-0.306206970039397,4.57373005245357,-4.38624938482081,0.706292855271531,0),
c("GO:0098840","protein transport along microtubule",0.00102330593773478,6.29560405345719,-0.781712183194048,2.39445168082622,-4.15116969303531,0.608816327864358,0.00606553),
c("GO:0042445","hormone metabolic process",0.11912192561817,-2.59932333570784,5.97312195483813,4.45869826840275,-3.61791651569513,0.956476101865797,0.03075045),
c("GO:0043651","linoleic acid metabolic process",0.00188918019274112,-1.91153268353234,-5.80089471466892,2.65991620006985,-3.25550860139878,0.987652127946439,0.03386545),
c("GO:0060336","negative regulation of interferon-gamma-mediated signaling pathway",0.000878303072063856,0.912812987218979,7.55195513862336,2.32837960343874,-3.10390072763946,0.894788741541073,0.12115344),
c("GO:0009251","glucan catabolic process",0.136928277522559,-4.34714429805278,-3.68138289535775,4.51919774408441,-3.05964538799872,0.985909990466915,0.12492643),
c("GO:0007218","neuropeptide signaling pathway",0.0819721914332804,-3.1723422527923,4.26449816000589,4.29637995377139,-3.29418273001038,0.939163593473231,0.16086415),
c("GO:0060013","righting reflex",0.00136716987632581,-4.69823141600724,2.65492574288902,2.51982799377572,-3.17830801158322,0.986179668377807,0.1872236),
c("GO:0015805","S-adenosyl-L-methionine transport",0.00176074908314688,6.21514101659149,3.17343102441649,2.62940959910272,-3.32360613294033,0.843654209392942,0.23883122),
c("GO:0015720","allantoin transport",0.000484723865242788,4.23780134560317,0.774558085796737,2.07188200730613,-3.18256173002017,0.845585828839705,0.24722283),
c("GO:0015855","pyrimidine nucleobase transport",0.0128762544715777,6.68179702148546,1.0346133691282,3.49262072204319,-3.38512477906643,0.797701641584638,0.28687488),
c("GO:0072334","UDP-galactose transmembrane transport",0.00306163193516599,5.42603053060789,1.6930999238571,2.86923171973098,-3.23376034926276,0.811785636565216,0.30312109),
c("GO:0033602","negative regulation of dopamine secretion",0.000290005731341839,0.173653337740208,7.3302508917582,1.85125834871908,-3.05611276737559,0.93827104356469,0.30888252),
c("GO:0015746","citrate transport",0.00647127074794218,5.3095697540338,2.28565028770444,3.19395897801919,-3.19119987391584,0.817848611740896,0.5354378),
c("GO:0099111","microtubule-based transport",0.058324295511863,6.221678606918,-2.45569186838678,4.14857180893216,-3.31472842369991,0.62792936419695,0.59563132),
c("GO:0060429","epithelium development",0.279511666806284,2.17898221223189,-6.87732030913693,4.82909783625819,-3.10354041739211,0.955684293389268,0.6050412),
c("GO:0099118","microtubule-based protein transport",0.00102330593773478,5.88553139633021,-0.815374942989983,2.39445168082622,-4.15116969303531,0.609563299604385,0.64774905),
c("GO:0015857","uracil transport",0.000658727304047892,7.15774886276768,1.47974800485151,2.20411998265592,-3.18256173002017,0.817571189789625,0.66828862));

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
ex <- one.data [ one.data$dispensability < 0.65, ];
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
ggsave("figures/DE_E1516_vs_E217_up_BP_REVIGO.png", width=8, height=8, dpi=1200);


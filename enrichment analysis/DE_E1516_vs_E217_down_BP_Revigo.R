# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800
# http://revigo.irb.hr/Results.aspx?jobid=1945704689
# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
library( ggrepel)

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0006695","cholesterol biosynthetic process",0.0148524363837213,-2.18873515028195,-5.22091778431574,3.55461028522616,-4.85514690693936,0.861745274040473,0),
c("GO:0030029","actin filament-based process",0.377363743500039,6.32712096637477,-3.78850754222606,4.95945639858514,-5.40307659172421,0.990226588339671,0),
c("GO:0080022","primary root development",0.000874160133044687,-4.69319534807439,4.179015370582,2.32633586092875,-4.6604444283948,0.842833436235332,0),
c("GO:0080036","regulation of cytokinin-activated signaling pathway",0.00311134720339602,-6.90440929285323,-0.769938108928796,2.87621784059164,-3.40243499217189,0.96450109400409,0),
c("GO:0071368","cellular response to cytokinin stimulus",0.0237100400067049,4.21185180546435,-5.54669052769081,3.75769962508774,-3.40243499217189,0.991743221456652,0.00782118),
c("GO:0072488","ammonium transmembrane transport",0.0708691148619071,0.519640741573832,7.78282282323375,4.23317385538094,-3.66342209228352,0.944076176791034,0.00848856),
c("GO:0010400","rhamnogalacturonan I \nside chain metabolic process",0.000227861646054302,5.53607477350318,0.782122146562825,1.7481880270062,-4.28040889770618,0.734553656600824,0.03190301),
c("GO:0090181","regulation of cholesterol metabolic process",0.00412636726309245,-5.87593765922476,-2.74155617439478,2.99869515831166,-3.35837029333796,0.964148620594852,0.1165392),
c("GO:0030032","lamellipodium assembly",0.00537753484688153,5.96230595973537,3.01471080821643,3.11360915107303,-4.17089285014785,0.803875878574667,0.21226024),
c("GO:0009263","deoxyribonucleotide biosynthetic process",0.184456073950467,-0.744616429563635,-5.67296420634953,4.64859417407893,-4.02465570701311,0.871206728857534,0.2830333),
c("GO:0030036","actin cytoskeleton organization",0.361612289349158,4.82810279901461,2.44984677644743,4.94093961624909,-4.12871633989144,0.741477196619817,0.31348757),
c("GO:0045026","plasma membrane fusion",0.0119316643752071,4.70273370307236,3.80793202796343,3.45954325828041,-3.11423146427531,0.759218322108861,0.33053434),
c("GO:0010075","regulation of meristem growth",0.00321906361789441,-4.50468851446755,3.16438132766063,2.89097959698969,-3.89705660804548,0.811789296630998,0.34057747),
c("GO:0010232","vascular transport",0.00309477544731934,-3.01891776155919,5.05988252771041,2.87390159786446,-3.23883597208328,0.802074520590397,0.34391436),
c("GO:0009395","phospholipid catabolic process",0.0516790213251157,-1.24970586677883,-5.0095077944348,4.09604055429543,-3.19200112479341,0.868896751040636,0.36943769),
c("GO:0009262","deoxyribonucleotide \nmetabolic process",0.260487290830259,-0.18566287496392,-5.35094321645584,4.79848490527716,-3.16385972731782,0.886206543946334,0.50439436),
c("GO:0097581","lamellipodium organization",0.00648784250401886,6.00499537476062,2.44704826664453,3.19506899646859,-4.17089285014785,0.802444929811298,0.57061658),
c("GO:0010233","phloem transport",0.00042672271897442,-3.4349650306672,5.04028971788467,2.01703333929878,-3.23883597208328,0.814634547227589,0.66243681),
c("GO:0031629","synaptic vesicle fusion to \npresynaptic active zone membrane",0.00159088858336095,2.83281484171003,3.84145399514596,2.5854607295085,-3.11423146427531,0.65693299651306,0.67378106),
c("GO:0099531","presynaptic process involved in \nchemical synaptic transmission",0.0349956031777321,-2.66889702022739,4.07099718083203,3.92675390518974,-3.00540896169249,0.780730773970911,0.68384721));

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
ggsave("figures/DE_E1516_vs_E217_down_BP_REVIGO.png", width=8, height=8, dpi=1200);

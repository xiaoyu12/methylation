# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800
#http://revigo.irb.hr/Results.aspx?jobid=1914903134

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
revigo.data <- rbind(c("GO:0015805","S-adenosyl-L-methionine transport",0.00176074908314688,3.13674137979025,4.17142491325915,2.62940959910272,-6.84340891257352,0.890080778422216,0),
c("GO:0042445","hormone metabolic process",0.11912192561817,0.0411180431766293,-5.89285402268519,4.45869826840275,-4.65797170826452,0.629991808468126,0),
c("GO:0048527","lateral root development",0.0052739613714023,3.35347986400274,-2.25995482774798,3.10516942799933,-3.47120160546549,1,0),
c("GO:0042365","water-soluble vitamin catabolic process",0.0294977258164842,-6.20698391513571,1.89197403111439,3.8525409857698,-4.02946141621083,0.793067718205775,0.02344119),
c("GO:0006081","cellular aldehyde metabolic process",0.427663166131772,6.38162493871825,-1.63712251720255,5.01379751314532,-3.44087091548831,0.956555175285919,0.04734933),
c("GO:0009767","photosynthetic electron transport chain",0.0301025949132829,4.39110380680688,-4.96507576933746,3.86135516019326,-3.15467993000543,0.963020378043479,0.0474202),
c("GO:0051224","negative regulation of protein transport",0.0229104527760053,-2.51697899603541,-6.55028850013132,3.74280365846917,-3.37252400265999,0.728731976405146,0.1472962),
c("GO:0007215","glutamate receptor signaling pathway",0.0495329789131861,-3.77285371852872,-5.04700024774429,4.07762222924452,-3.18879113334673,0.877833651183322,0.15521446),
c("GO:0006839","mitochondrial transport",0.136402124267124,5.10380866123969,2.44967832106416,4.51752578363382,-4.20565946895543,0.923287225901378,0.16031369),
c("GO:0090487","secondary metabolite catabolic process",2.07146950958457E-05,-6.54129018232215,-0.734392844339635,0.778151250383644,-3.2039582459174,0.896857408261797,0.22354787),
c("GO:1990818","L-arginine transmembrane export from vacuole",0.000219575768015964,2.62729544511743,5.73090060203222,1.73239375982297,-4.0574031964141,0.787422464886496,0.26897437),
c("GO:0006598","polyamine catabolic process",0.00473537929891032,-4.06273991031246,2.80922474072456,3.05842602445701,-3.92662416341722,0.766218733967936,0.29068678),
c("GO:1901264","carbohydrate derivative transport",0.246695446835445,4.09988573439285,4.23630873193166,4.77485988645246,-3.62886703268881,0.877303831320046,0.31610943),
c("GO:0072526","pyridine-containing compound catabolic process",0.0250274946148007,-5.45234685194122,2.70609882251959,3.78118072093726,-3.20300858291416,0.786998504408148,0.39126925),
c("GO:0009407","toxin catabolic process",0.00372450217823305,-5.55400305924304,-0.612608318706118,2.95424250943932,-3.2039582459174,0.847061224577207,0.44205197),
c("GO:0006843","mitochondrial citrate transmembrane transport",0.000679441999143737,3.22349791207222,5.55825245768291,2.21748394421391,-3.47104208168567,0.842193568321882,0.44976652),
c("GO:0009111","vitamin catabolic process",0.0324972136663627,-5.79813528513132,1.63834744348044,3.89459294792296,-3.64114688683226,0.812603785108215,0.52792371),
c("GO:0010817","regulation of hormone levels",0.198823786468946,-0.411706397257749,-5.72380754329917,4.6811688489294,-3.12471650450369,0.756406648316712,0.53273337),
c("GO:0006595","polyamine metabolic process",0.126662074633058,-2.85342495411135,6.77375860713519,4.48535226124114,-3.07812471782251,0.825833025068469,0.64303899),
c("GO:0042822","pyridoxal phosphate metabolic process",0.0940240010400434,-5.77296461850766,4.4701056830497,4.3559493227878,-3.02747635833075,0.795289361974714,0.68991365));

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
ggsave("figures/DE_shared_up_top_BP_REVIGO.png", width=8, height=8, dpi=1200);

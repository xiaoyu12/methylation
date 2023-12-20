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
revigo.data <- rbind(c("GO:0009378","four-way junction helicase activity",0.0630908438695814,-1.42994760484879,6.83221185852912,4.26147698862137,-3.93330050019013,0.768977994530058,0),
c("GO:0032135","DNA insertion or deletion binding",0.00232210795707957,6.04987426767515,0.996564809996078,2.82801506422398,-2.36491487417872,0.715853519406694,0),
c("GO:0009702","L-arabinokinase activity",0.000310996601394585,-5.42178412843178,1.29555216708676,1.95904139232109,-2.4325145879317,0.917989766510501,0.01652816),
c("GO:0004610","phosphoacetylglucosamine mutase activity",0.005034689424799,2.4822170613288,-5.91485828873707,3.16375752398196,-2.38684735990483,0.886955203311552,0.01919367),
c("GO:0000257","nitrilase activity",0.0134074090378999,-3.44324046690925,-4.42389202448193,3.58894364274001,-2.4084035768878,0.915958547165649,0.02034804),
c("GO:0004349","glutamate 5-kinase activity",0.0297001754331828,-5.22896392525123,2.5707392930743,3.93429640681941,-2.38375963687872,0.904962892368856,0.25239298),
c("GO:0003896","DNA primase activity",0.0488713881502624,-3.19851841824119,4.76708478882808,4.150572247669,-2.14403950916006,0.775232499276906,0.338961),
c("GO:0016815","hydrolase activity, \nacting on carbon-nitrogen (but not peptide) bonds, \nin nitriles",0.0135318076784577,-4.03615784213081,-4.02755224224233,3.59295357154787,-2.4084035768878,0.915947423378271,0.44953121),
c("GO:0001006","RNA polymerase III type 3 \npromoter sequence-specific DNA binding",0.00193854548202625,5.33833031582121,2.37262072878405,2.74973631556906,-2.33290291426046,0.721648522033571,0.46038452),
c("GO:0034458","3'-5' RNA helicase activity",0.00221844242328137,-0.79765924369883,7.03950145387667,2.80821097292422,-2.04601570210968,0.802379893559345,0.47511493),
c("GO:0001016","RNA polymerase III transcription regulatory \nregion sequence-specific DNA binding",0.0022529976012141,6.01678088307557,2.28114628540267,2.81491318127507,-2.33290291426046,0.720118451139543,0.51687455),
c("GO:0032138","single base insertion or deletion binding",0.00158262714931911,5.48121370872687,0.89778974515323,2.66181268553726,-2.36491487417872,0.719516107326156,0.62153836),
c("GO:0004614","phosphoglucomutase activity",0.0138151601375061,1.96780202886027,-5.94596247709284,3.60195140413352,-2.36081629071656,0.88581000910146,0.67067219),
c("GO:0000976","transcription cis-regulatory region binding",0.373271943064954,5.64814085894424,1.71781712588238,5.03351623427969,-2.08952990710244,0.667056553153636,0.68367933));

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
ex <- one.data [ one.data$dispensability < 0.35, ];
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
ggsave("figures/DMPs_MF_REVIGO.png", width=8, height=8, dpi=1200);

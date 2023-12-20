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
revigo.data <- rbind(c("GO:0008519","ammonium transmembrane transporter activity",0.0591031763361442,4.78108498250751,4.44621270850827,4.23312307852108,-3.34083763711449,0.876815049595946,0),
c("GO:0051015","actin filament binding",0.133376075784758,-5.51847228681141,-2.32831323137556,4.58657605337363,-6.97552277184929,1,0),
c("GO:0061731","ribonucleoside-diphosphate reductase activity",0.0684641740381212,-4.4027507185081,4.27240031886151,4.29697215865653,-4.59579901465999,0.577826481434966,0),
c("GO:0004623","phospholipase A2 activity",0.038152371955529,5.03179178193271,-3.61393974746887,4.04304774280412,-3.50790873121506,0.805935553203649,0.02186161),
c("GO:0004611","phosphoenolpyruvate carboxykinase activity",0.0518224003457176,-0.542584593050875,-5.31895938041014,4.17603334926402,-3.2515410609627,0.970895000737993,0.02231239),
c("GO:0016682","oxidoreductase activity, \nacting on diphenols and related substances as donors, \noxygen as acceptor",0.0418186763341918,-2.94458574592852,6.04992295737433,4.08289303328818,-3.39079843509358,0.834594217949343,0.23694347),
c("GO:0005216","ion channel activity",1.06305549392256,5.16632588276712,3.59907737925633,5.4880442143777,-3.15223470245304,0.876815049595946,0.57254374));

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
ggsave("figures/DE_E1516_vs_E217_down_MF_REVIGO.png", width=8, height=8, dpi=1200);

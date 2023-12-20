

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0015805","S-adenosyl-L-methionine transport", 0.002,-4.295,-4.612, 2.326,-6.8434,0.493,0.000),
c("GO:0042445","hormone metabolic process", 0.090, 1.697, 5.690, 4.064,-4.6580,0.823,0.000),
c("GO:0042365","water-soluble vitamin catabolic process", 0.029, 5.213,-3.674, 3.570,-4.0295,0.620,0.043),
c("GO:0006839","mitochondrial transport", 0.182,-4.483, 1.678, 4.369,-4.2057,0.635,0.165),
c("GO:0089707","L-lysine transmembrane export from vacuole", 0.002,-3.632,-2.111, 2.322,-4.0574,0.298,0.259),
c("GO:0051181","cofactor transport", 0.109,-5.262,-1.678, 4.146,-5.1332,0.506,0.303),
c("GO:0006598","polyamine catabolic process", 0.002, 5.486,-1.658, 2.405,-3.9266,0.722,0.304),
c("GO:0051182","coenzyme transport", 0.018,-4.742,-3.794, 3.359,-5.7656,0.475,0.693),
c("GO:0042820","vitamin B6 catabolic process", 0.000, 4.994,-4.264, 0.845,-4.0295,0.626,0.711));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );

#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 1, ]; 
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
s <- 5
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/s,max(one.data$plot_X)+one.x_range/s);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/s,max(one.data$plot_Y)+one.y_range/s);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

ggsave("DE_shared_up_top_GO_REVIGO.png", device="png", dpi=600);

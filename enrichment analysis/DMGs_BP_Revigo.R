# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800
# http://revigo.irb.hr/Results.aspx?jobid=1947760782

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
revigo.data <- rbind(c("GO:0051354","negative regulation of oxidoreductase activity",0.00354635580040878,1.16833871211128,-6.20387792092927,2.9329808219232,-3.19233965229192,0.808604335202994,0),
c("GO:0060712","spongiotrophoblast layer development",0.00211704183879543,-4.82604364675679,-2.60726774359003,2.70926996097583,-2.07605265261839,0.949882851427776,0),
c("GO:0071277","cellular response to calcium ion",0.0157473112118619,-5.66989212277966,2.60840939522009,3.58001211252942,-2.50878722363399,0.995299359645842,0),
c("GO:1904715","negative regulation of \nchaperone-mediated autophagy",0.000310720426437685,5.37106059761772,4.02022552124552,1.88081359228079,-2.85480855124118,0.778626862273426,0.09709413),
c("GO:0055088","lipid homeostasis",0.0508504335212819,5.6052108832275,-3.42000036513529,4.08902150079501,-2.44388606641999,0.764153182439323,0.12485255),
c("GO:0048831","regulation of shoot system development",0.0134686947513188,1.50646089540394,-2.48975237289947,3.51215053692203,-2.00279454010732,0.864548824856279,0.13495216),
c("GO:0034250","positive regulation of cellular amide \nmetabolic process",0.0915175229334461,6.243102767646,0.436194643584057,4.3442153756565,-2.1267972163376,0.830150796239848,0.18685739),
c("GO:0031507","heterochromatin assembly",0.0349374047486533,2.23012111322964,3.74448263202255,3.92603359667884,-2.63774040385431,0.478887663212671,0.38970448),
c("GO:0044782","cilium organization",0.171923683417481,-1.39017778963761,6.22231847535363,4.61803763165873,-2.17715662335989,0.756185347431888,0.42731391),
c("GO:0002832","negative regulation of response to biotic stimulus",0.0226038752885868,3.46331866827891,1.66791263607526,3.73695395378315,-2.2184266145008,0.808629530237222,0.43714264),
c("GO:0031017","exocrine pancreas development",0.00429208482385922,-4.47868264200269,-3.39956140338982,3.01577875638904,-2.06486451013155,0.949882851427776,0.45685194),
c("GO:1904714","regulation of \nchaperone-mediated autophagy",0.000716728450316259,6.76765701151176,4.0444516926901,2.2405492482826,-2.23107358990402,0.818427895680583,0.58648639),
c("GO:0042632","cholesterol homeostasis",0.0184070780621684,5.6546180776144,-4.10119627862732,3.64777405026883,-2.58688500832399,0.744206326816836,0.61728692),
c("GO:0071103","DNA conformation change",1.12146873485497,-0.767182934778698,6.49981606485207,5.432480234003,-2.1651518317695,0.677213061832958,0.6251554));

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
ex <- one.data [ one.data$dispensability < 0.75, ];
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
ggsave("figures/DMGs_BP_REVIGO.png", width=8, height=8, dpi=1200);

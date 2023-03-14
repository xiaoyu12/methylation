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
revigo.data <- rbind(c("GO:0008569","minus-end-directed microtubule motor activity",0.0398904974055454,-1.80722765931445,-1.80214634517111,4.0623939372532,-3.68973668293612,1,0),
c("GO:0008810","cellulase activity",0.0479107542037324,-5.49112939038624,-2.99540313486925,4.14195119586275,-6.20228346785532,0.973182505399619,0),
c("GO:0015391","nucleobase:cation symporter activity",0.000749847361140277,5.20587173431941,4.05394028419751,2.3384564936046,-3.18256173002017,0.662783348023958,0),
c("GO:0045504","dynein heavy chain binding",0.00357646091603772,-4.09812338279407,4.52760002602165,3.01535975540921,-5.6013625727694,0.897116460465566,0),
c("GO:0015020","glucuronosyltransferase activity",0.0428311430476209,0.541208963513703,-5.79914153402321,4.09328156756725,-3.05018265886471,0.92108952886723,0.02151906),
c("GO:0016594","glycine binding",0.00506233356714519,6.65342057874868,-3.57440897471021,3.16613397030511,-3.41636664357001,0.987396846869623,0.02700489),
c("GO:0016435","rRNA (guanine) methyltransferase activity",0.0656168273764641,3.21934568063216,-5.29032493506088,4.27852496473702,-3.0117336151904,0.958528699310244,0.15890671),
c("GO:0008484","sulfuric ester hydrolase activity",0.165225583285356,-6.21252764326452,-1.12868707592472,4.67957324282562,-5.16163595103143,0.972642339468988,0.18307853),
c("GO:0045503","dynein light chain binding",0.00263656007626742,-2.3537210154301,6.24377907045403,2.88309335857569,-3.32562323804185,0.898257040478614,0.28369003),
c("GO:0016934","extracellularly glycine-gated chloride channel activity",0.00336567433064806,5.92294945704962,2.0170942352456,2.98900461569854,-3.17830801158322,0.74519345812657,0.31030329),
c("GO:0045505","dynein intermediate chain binding",0.0300077165167842,-3.46895925148509,5.68416599127521,3.93876982278312,-4.46388987587942,0.889025668165549,0.32128826),
c("GO:0022857","transmembrane transporter activity",9.00077033858165,5.18360214117796,2.66273520029867,6.41576671362888,-3.18256173002017,0.763856430764102,0.34258057),
c("GO:0051959","dynein light intermediate chain binding",0.0264519887075061,-2.97850736433822,4.97007256174177,3.88400192476879,-3.4300192448464,0.889499854966985,0.36061064),
c("GO:0080116","glucuronoxylan glucuronosyltransferase activity",0.000297174530221492,-0.360180580193669,-5.92226319089949,1.93951925261862,-3.30232683298143,0.926505615976374,0.52096031),
c("GO:0005274","allantoin:proton symporter activity",0.00040429558181296,4.5835409121323,4.49296900238839,2.07188200730613,-3.18256173002017,0.743795902364129,0.53725068),
c("GO:0005350","pyrimidine nucleobase transmembrane transporter activity",0.0106533613566612,5.78396570887083,3.86373850504444,3.48911436937892,-3.18256173002017,0.652097299534101,0.6840284));

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
ggsave("figures/DE_E1516_vs_E217_up_MF_REVIGO.png", width=8, height=8, dpi=1200);

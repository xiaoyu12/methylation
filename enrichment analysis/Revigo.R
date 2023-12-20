# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
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
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000095","S-adenosyl-L-methionine transmembrane transporter activity",0.00146513954434782,-4.22552943106912,4.18391338974828,2.62838893005031,-4.89268074933159,0.944440188965492,0),
c("GO:0015036","disulfide oxidoreductase activity",0.325506320608539,-4.194785427027,-4.02818843992382,4.97405090279288,-5.08880611969258,0.644270057074587,0),
c("GO:0070402","NADPH binding",0.0323747462051763,1.70156851874585,-6.80634939052463,3.97173959088778,-3.28581133117432,0.993758992895778,0),
c("GO:0016277","[myelin basic protein]-arginine N-methyltransferase activity",0.000252252798908941,5.69011688377037,-0.617307616916006,1.86923171973098,-3.99618447893561,0.741459723394491,0.01779896),
c("GO:0008157","protein phosphatase 1 binding",0.00863188344759636,3.69791975512997,5.79403204720422,3.39776625612645,-3.00543934118677,0.993758992895778,0.03125388),
c("GO:0072349","modified amino acid transmembrane transporter activity",0.0133797648955537,-2.64097327463587,5.07522698304847,3.58804749698608,-3.33020737934869,0.944132542045571,0.26726342),
c("GO:0016667","oxidoreductase activity, acting on a sulfur group of donors",0.574939416998169,-5.19521939947339,-2.94723276224065,5.22111156086213,-8.95273909427272,0.816149658316953,0.32792302),
c("GO:0035241","protein-arginine omega-N monomethyltransferase activity",0.00202493342685807,5.58210217090679,-1.31793951303418,2.76863810124761,-3.44501669650421,0.716519788817192,0.6527751));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
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


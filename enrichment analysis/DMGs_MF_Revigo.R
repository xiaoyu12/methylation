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
revigo.data <- rbind(c("GO:0001786","phosphatidylserine binding",0.00880465933726002,-5.81632941078689,-2.86359707269451,3.40636983546927,-4.16399439370572,0.817091045115801,0),
c("GO:0004656","procollagen-proline 4-dioxygenase activity",0.0108572369064643,-2.56163646527288,5.98928546230546,3.49734438101758,-2.4497297226091,0.812187320908162,0),
c("GO:0061134","peptidase regulator activity",0.218115738629195,5.48570136507636,-3.53635316624251,4.8001807508985,-2.29762048435824,0.92847015876258,0),
c("GO:0016856","racemase and epimerase activity, \nacting on hydroxy acids and derivatives",0.0139050036001312,3.27067184783555,4.62124983246448,3.60476588470389,-2.26838709846644,0.921315248410953,0.01840088),
c("GO:0016409","palmitoyltransferase activity",0.080226756606423,-2.81520932645886,-7.20581588113295,4.36582480685936,-2.21470775133807,0.982273692422163,0.02069977),
c("GO:0051880","G-quadruplex DNA binding",0.00255017213143559,-5.9631271307723,3.18965733121455,2.86864443839483,-2.50731511268998,0.959256622302953,0.02729677),
c("GO:0019905","syntaxin binding",0.0254256999229039,0.773729277541532,-7.67749652043218,3.86681880292605,-2.02469622396295,0.979293013139269,0.03086034),
c("GO:0140666","annealing activity",0.0110887565986136,5.84896429381661,1.50730520753869,3.50650503240487,-2.32192504672345,0.957143875588585,0.1355803),
c("GO:0015278","calcium-release channel activity",0.0111129452231665,3.05575811008032,-3.48717770785629,3.50745106090197,-2.00490783539893,0.940628566788782,0.1566335),
c("GO:0052596","phenethylamine:\noxygen oxidoreductase (deaminating) activity",0.00639270791755535,-2.11258603936348,4.23165034308475,3.2674064187529,-2.18905263551172,0.734847259200927,0.19002786),
c("GO:0010521","telomerase inhibitor activity",0.00158953818490566,4.46722536040956,-3.70023760428496,2.66370092538965,-2.03600503698869,0.905359655916413,0.49907835),
c("GO:0032934","sterol binding",0.0183522550000738,-5.39574095269438,-2.32201395829522,3.72525806635996,-2.42429528820131,0.770030126945201,0.54847516),
c("GO:0005496","steroid binding",0.0293788122784084,-6.07006590290591,-2.17792652894657,3.92957217907655,-2.0004824103384,0.802150444115648,0.58796581),
c("GO:0008692","3-hydroxybutyryl-CoA epimerase activity",0.00513489944080392,2.93913414691007,5.02188385869009,3.17231096852195,-2.26838709846644,0.922042249462346,0.59340094));

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
ggsave("figures/DMGs_MF_REVIGO.png", width=8, height=8, dpi=1200);

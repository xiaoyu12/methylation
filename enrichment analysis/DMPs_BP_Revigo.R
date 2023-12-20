# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800
# http://revigo.irb.hr/Results.aspx?jobid=1947977764
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
revigo.data <- rbind(c("GO:0007585","respiratory gaseous exchange by respiratory system",0.00710928335689423,3.82531907644155,-5.03702949969074,3.23477029516092,-2.48558155136539,0.913097265491011,0),
c("GO:0061744","motor behavior",0.00223304413133216,6.53524478496151,0.688282947887486,2.73239375982297,-2.44906346271243,0.966592455840543,0),
c("GO:0070129","regulation of mitochondrial translation",0.00544382187118824,-2.10963130624591,3.10702837671239,3.11892575282578,-2.7180579621706,0.885480877624097,0),
c("GO:0071215","cellular response to abscisic acid stimulus",0.0320373474352349,-5.11511594079141,-4.43383472213614,3.88840416773705,-3.38392804915916,0.900526780916418,0),
c("GO:0080186","developmental vegetative growth",0.000774729596584627,5.51517480728556,-2.64547720282899,2.27415784926368,-3.08361973055843,0.883943593223802,0),
c("GO:0019255","glucose 1-phosphate metabolic process",0.000261005158207655,0.884591881331911,-7.63813165087057,1.80617997398389,-2.38684735990483,0.990321682712289,0.00521111),
c("GO:0008298","intracellular mRNA localization",0.00780115417309547,-0.473286816485504,6.99183731679516,3.27508089845686,-2.48270755218371,0.904517935463228,0.00622254),
c("GO:0032981","mitochondrial respiratory chain complex I assembly",0.0392709189627042,4.47987521976518,4.55312806783694,3.97680833733807,-2.31855278293552,0.898415818210584,0.00685551),
c("GO:0042816","vitamin B6 metabolic process",0.131003874725147,-0.921122908369027,-6.48763138661062,4.49998933433418,-2.02493826981282,0.969096690495953,0.03299276),
c("GO:1904058","positive regulation of sensory perception of pain",0.00033143512153353,-6.43179951436244,2.18060473501675,1.90848501887865,-2.67473785223923,0.764963999244552,0.10634632),
c("GO:0048872","homeostasis of number of cells",0.0519565982394001,-0.692262909092514,0.597834335627676,4.09836679643933,-2.28217796462103,0.896071945792099,0.12808561),
c("GO:0019985","translesion synthesis",0.0281222700621201,-3.37007179310568,-4.57459638447563,3.83180580867439,-2.36623494777214,0.893646651630613,0.24109694),
c("GO:0045945","positive regulation of transcription by RNA polymerase III",0.00331020827631614,-4.73387117708936,3.81972014933691,2.90308998699194,-2.12028899498221,0.850531606896067,0.32949063),
c("GO:1905273","positive regulation of proton-transporting \nATP synthase activity, rotational mechanism",0.00047643798720445,-5.84686647544338,4.12407336159151,2.06445798922692,-2.0137303535048,0.830190951441745,0.34329911),
c("GO:0003016","respiratory system process",0.00367478691000302,3.35282442275037,-4.63086468617346,2.9484129657786,-2.44906346271243,0.915421245278798,0.36069144),
c("GO:0070417","cellular response to cold",0.010386348121057,-4.32522881638233,-5.28236354533547,3.39932753215868,-2.29275896212231,0.925417982151836,0.38920987),
c("GO:0021510","spinal cord development",0.022168866691574,4.53038623494809,-4.0054958605115,3.72851610475977,-2.27732086046529,0.852019965416115,0.39677299),
c("GO:0006403","RNA localization",0.201839846074901,-0.890587514979189,7.22496740387355,4.68770727962482,-2.03613008741456,0.934147303324885,0.41091697),
c("GO:0000280","nuclear division",0.233019605133168,4.06616187871387,5.11627423741512,4.75009164252212,-2.29986296731142,0.903643315078259,0.44025229),
c("GO:0009952","anterior/posterior pattern specification",0.0542807870291539,4.24301583113565,-3.78595959370749,4.11737074102091,-2.23056849881708,0.852329426925788,0.51404865),
c("GO:0045579","positive regulation of B cell differentiation",0.00286277086224587,-5.82252998185445,2.03436919465077,2.84010609445676,-2.62754029128666,0.704322564974176,0.51937863),
c("GO:0014876","response to injury involved in \nregulation of muscle adaptation",0.000153288743709258,-6.14943156998507,-1.00880799414881,1.57978359661681,-2.44906346271243,0.757362290077082,0.53364426),
c("GO:0038179","neurotrophin signaling pathway",0.00536096309080485,-4.42712010801212,-2.83629075681047,3.11226976841727,-2.04790641286915,0.812467423610988,0.54506726),
c("GO:0032508","DNA duplex unwinding",0.658317153084994,3.94521014194227,4.73245827900956,5.20112936343455,-2.1662906083894,0.898651232364464,0.55221942),
c("GO:0010257","NADH dehydrogenase complex assembly",0.0401657937908447,4.45260991658905,4.12260514499979,3.98659260682221,-2.31855278293552,0.908704579186803,0.56221838),
c("GO:0036343","psychomotor behavior",0.000609012035817862,6.77477116440575,0.503766793737788,2.17026171539496,-2.44906346271243,0.966592455840543,0.58959811),
c("GO:0071033","nuclear retention of pre-mRNA at the site of transcription",0.00360021400765797,-1.56212783336688,5.17967011379862,2.93951925261862,-2.69587631993226,0.780999203175289,0.62902403),
c("GO:0038180","nerve growth factor signaling pathway",0.00102330593773478,-4.85929887029073,-3.39019164057942,2.39445168082622,-2.04790641286915,0.821564396893206,0.63797361),
c("GO:1903706","regulation of hemopoiesis",0.0584982989506681,-5.85595414367121,1.01025568815437,4.14986545302626,-2.05331265942471,0.731136825156186,0.65884192));

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
ggsave("figures/DMPs_BP_REVIGO.png", width=8, height=8, dpi=1200);

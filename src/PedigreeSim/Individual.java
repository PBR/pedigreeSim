/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Individual extends class Genotype:
 * an Individual has a name and can have parents.
 * Individuals can do meiosis tot produce gametes and they can mate to
 * produce offspring.
 * @author Roeland Voorrips
 */

public class Individual extends Genotype {
    private String indivName;
    private Individual[] parents;

    /**
     * parconfig stores the meiotic configurations that produced the two
     * Gametes from which this Individual originated.
     * The first index is for the parent (0 or 1), the second for the chromosome
     */
    ChromConfig[][] parconfig;

    /* The next 3 data structures are created and calculated by calcGenomehs().
     * They will contain information on the genome numbers (see below) at the
     * two ends (heads and tails) of the haplostructs. They will be created
     * the first time they are needed (when doMeiosis is called).
     * These data structures will only contain valid information
     * if ploidy>2 and only for chromosomes whose prefPairingProb > 0.0.
     */
    
    /**
     * There are in principle ploidy/2 different diploid genomes in a polyploid
     * population, and a founder allele fa belongs to genome fa % (ploidy/2).
     * genomenr[chrom][side][h] will contain the genome number corresponding
     * to the founderallele at side 0 or 1 (head or tail) of haplostruct[h]
     * of chromosome chrom.
     */
    int[][][] genomenr; 
    
    /**
     * genomehs will contain (only for polyploids, and only for chromosomes
     * where prefPairingProb>0), for each "genome number" (see above), which 
     * haplostruct ends have this number. 
     * genomehs has four index levels, where
     * genomehs.get(chrom).get(side).get(genomenr) returns an ArrayList<Integer>
     * that contains the indices to all haplostructs of chromosome chrom
     * that at side 0 or 1 (head or tail) have a founderallele belonging to
     * that genomenr.
     */
    ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>> genomehs;
    
    /**
     * possiblePairCount will contain the maximum number of preferential
     * pairs at both sides of each chromosome. A preferential pair
     * is a pair where both members have the same genome number (see above).
     * If each genomenr occurs exactly twice, or generally if no genome
     * occurs an odd number of times, there are ploidy/2 preferential pairs 
     * possible, else less preferential pairs are possible.
     */
    int[][] possiblePairCount;


    /**
     * Constructor to create a basic individual with only a name;
     * without parents, so it is a founder, and without haplostruct
     * @param name
     * @param popdata
     * @throws Exception
     */
    public Individual(String name, PopulationData popdata) throws Exception  {
        super(popdata);
        this.tools = popdata.tools;
        this.rand = popdata.tools.rand;
        this.indivName = name;
        this.parents = new Individual[] {null, null};
    }

    /**
     * Constructor to create a founder: its [chromCount][ploidy] HaploStructs
     * each have only one segment, with the founderalleles  =
     * (lastfounderallele+1) .. (lastfounderallele+popdata.ploidy)
     * and its parents are set to {null,null}
     * @param name
     * @param lastfounderallele integer, the last founder allele used up to now
     * @param popdata
     * @throws Exception
     */
    public Individual(String name, int lastfounderallele,
            PopulationData popdata) throws Exception  {
        this(name, popdata);
        this.setHaploStruct(Individual.calcFounderHaplostruct(lastfounderallele, popdata));
    }

    /**
     * Constructor with parents explicitly specifying the haplostruct
     * @param name
     * @param haplostruct
     * @param parents
     * @param popdata
     * @throws Exception
     */
    public Individual(String name, HaploStruct[][] haplostruct,
            Individual[] parents, PopulationData popdata)
            throws Exception {
        this(name, popdata);
        setParents(parents);
        setHaploStruct(haplostruct);
    }

    /**
     * Constructor specifying the parents but not the haplostruct.
     * This would be the normal constructor for generating non-founders in
     * a pedigree before their haplotypes have been calculated.
     * @param name
     * @param parents
     * @param popdata
     * @throws Exception
     */
     public Individual(String name, Individual[] parents, PopulationData popdata)
            throws Exception {
        super(popdata);
        this.indivName = name;
        setParents(parents);
        this.haplostruct = null;
        //a call to setHaploStruct must be made before the individual
        //is used to generate gametes
    }

    protected final void setHaploStruct(HaploStruct[][] haplostruct) throws Exception {
        if (checkHaploStruct(popdata.ploidy, haplostruct)) {
            this.haplostruct = haplostruct;
        } else {
            throw new Exception("Individual setHaploStruct: invalid haplostruct");
        }
    }

    /**
     * calcFounderHaplostruct generates the haplostructs for a new founder
     * (and in the process assigns (ploidy) new founderalleles)
     * @param lastfounderallele the last founder allele assigned up to now
     * @param popdata
     * @return
     */
    public static HaploStruct[][] calcFounderHaplostruct(int lastfounderallele,
            PopulationData popdata) {
        HaploStruct[][] haplostruct;
        haplostruct = new HaploStruct[popdata.chromCount()][popdata.ploidy];
        for (int p=0; p<popdata.ploidy; p++) {
            lastfounderallele++;
            for (int c=0; c<popdata.chromCount(); c++) {
                try {
                    haplostruct[c][p] = new HaploStruct(popdata.getChrom(c),
                            lastfounderallele);
                } catch (Exception ex) {
                    Logger.getLogger(Individual.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        return(haplostruct);
    } //calcFounderHaplostruct

    public String getIndivName() {
        return indivName;
    }

    public void setIndivName(String indivName) {
        this.indivName = indivName;
    }

    public Individual[] getParents() {
        return parents;
    }

    /**
     * setParents is used to assign parents (only used in the Individual
     * constructors); it checks whether there are two parents, both non-null
     * @param parents
     * @throws Exception 
     */
    private void setParents(Individual[] parents) throws Exception {
        if (parents.length != 2) {
            throw new Exception("Individual.setParents: parents.length!=2");
        }
        if (parents[0]==null || parents[1]==null) {
            throw new Exception("Individual.setParents: parents may not be null");
        }
        this.parents = parents;
    }

    /**
     * 
     * @return true if this individual has no parents, else false
     */
    public boolean isFounder() {
        return parents[0]==null;
    }


    /**
     * mating obtains two Gametes,
     * the first derived from this parent, the second from otherParent;
     * then thisGamete.fertilization(offspringName, otherGamete) is called.
     * In the offspring the haplostruct from thisGamete are added first
     * and those from otherGamete last; the order of the haplostruct are
     * not changed within or between the gametes.
     * @param otherParent
     * @param zygote the (already initialized) Individual produced by this
     *     mating; its parents must already be specified
     * @throws Exception
     */
    public void mating(Individual otherParent, Individual zygote)
            throws Exception {
        if (otherParent==null ||
                !otherParent.getPopdata().equals(popdata)) {
            throw new Exception("Other parent invalid in mating");
        }
        if (zygote==null ||
            !zygote.getPopdata().equals(popdata))
            throw new Exception("zygote is invalid in mating");
        Individual[] zyparents = zygote.getParents();
        if (zyparents == null) {
            zyparents = new Individual[2];
            zyparents[0] = this;
            zyparents[1] = otherParent;
            zygote.setParents(zyparents);
        }
        else if (!zygote.parents[0].equals(this) ||
                 !zygote.parents[1].equals(otherParent))
            throw new Exception("zygote has invalid parents in mating");
        zygote.conception();
    }

    /**
     * conception creates the genotype (haplostruct) of this Individual from
     * the fusion of the two parental Gametes.
     * The parents of this Individual must have been specified earlier.
     * The haplostruct from the Gamete from parent[0] are added first
     * and those from parent[1] last; the order of the haplostruct are
     * not changed within or between the gametes.
     * The configurations of the parental meioses are also stored in
     * this Individual.
     * @throws Exception
     */
    public void conception() throws Exception {
        if (parents==null || parents[0]==null || parents[1]==null)
            throw new Exception("parents are null in conception");
        Gamete par0gamete = parents[0].doMeiosis().get(0);
        Gamete par1gamete = parents[1].doMeiosis().get(0);
        par0gamete.fertilization(par1gamete, this);
    }

    /**
     * doMeiosis produces four gametes.
     * If ploidy==2, each gamete contains one HaploStruct per chromosome;
     * if ploidy==4, each gamete contains two HaploStruct per chromosome, etc.
     * @return an ArrayList of size 4: 4 Gametes.
     *   each element is a HaploStruct[][], with
     *   first index is chromosome number, and
     *   second index is number of HaploStruct (0 for haploid gamete,
     *   0 and 1 for diploid gamete, etc.)
     *   The gametes are ordered as a tetrad, with the first two deriving
     *   from the "top" side of the first meiotic division and the last two
     *   from the "bottom" side. For each chromosome the decision which
     *   centromeres go to the top and which to the bottom side is
     *   random.
     *   So although the order of the gametes is meaningful, picking one
     *   Gamete (e.g. the first) yields a completely random gamete
     * @throws Exception
     */
    public ArrayList<Gamete> doMeiosis() throws Exception {
        ArrayList<Gamete> gametes = new ArrayList<>();
        ChromConfig[] allchromconf = new ChromConfig[popdata.chromCount()];
        //Gamete gam;
        HaploStruct[][] hs = new HaploStruct[popdata.chromCount()][2*popdata.ploidy];
        //hs[c] will contain the 2n haplostructs to go into the 4 alleles
        //in the correct sequence (first n/2 for gamete 0 etc)
        HaploStruct[] origChrom; //copy of haplostruct[c] in order of chrconf.chromseq
        HaploStruct[] qh; // = new HaploStruct[8]; //result of each quadrivalent
        HaploStruct[] dh; // = new HaploStruct[4]; //result of each bivalent
        int p2 = popdata.ploidy/2;
        for (int c=0; c<popdata.chromCount(); c++) {
            calcChromConfig(c);
            allchromconf[c] = new ChromConfig(chrconf);
            //first do the quadrivalents:
            origChrom = new HaploStruct[4];
            for (int v=0; v<chrconf.quad; v++) {

                for (int h=0; h<4; h++) {
                    origChrom[h] = haplostruct[c][chrconf.chromseq[4*v+h]];
                }    
                qh = new Quadrivalent(origChrom).doMeiosis();
                for (int g=0; g<4; g++) {
                    hs[c][2*v+p2*g] = qh[2*g];
                    hs[c][2*v+p2*g+1] = qh[2*g+1];
                }
            }
            //then do the bivalents:
            int firstbi = chrconf.quad*4;
            if (firstbi < popdata.ploidy) {
                origChrom = new HaploStruct[2];
                for (int v=0; v<(popdata.ploidy-firstbi)/2; v++) {
                    for (int h=0; h<2; h++) {
                        origChrom[h] = haplostruct[c]
                                [chrconf.chromseq[firstbi+2*v+h]];
                    }
                    dh = new Bivalent(origChrom).doMeiosis();
                    for (int g=0; g<4; g++) {
                        hs[c][firstbi/2+v+p2*g] = dh[g];
                    }
                }
            } 
            //randomize hs[c] within each gamete:
            //System.out.println("randomize hs[c] within each gamete for chrom "+c);
            HaploStruct[] shufhs = new HaploStruct[popdata.ploidy/2];
            for (int g=0; g<4; g++) {
                int[] shuf = tools.shuffleArray(Tools.seq(g*popdata.ploidy/2,
                        ((g+1)*popdata.ploidy/2)-1));
                //make a shuffled copy of the hs for this gamete in shufhs
                for (int h=0; h<popdata.ploidy/2; h++) {
                    shufhs[h] = hs[c][shuf[h]];
                }
                //copy the shuffled version back into hs:
                System.arraycopy(shufhs, 0, hs[c], g*popdata.ploidy/2, popdata.ploidy/2);
            }
        } //for c
        //finally put hs into gametes:
        //now we have the four gametes represented in hs as 4*(ploidy/2) HaploStructs
        //per chromosome; these are now converted into Gametes:
        
        //Modification 20200904 according to Soichiro Nishiyama,
        //to simulate unreduced (FDR and SDR) gametes, (test mode only!)
        if (popdata.unreducedGametes == 0) {
            //normal: 4 haploid gametes
            for (int g=0; g<4; g++) {
                gametes.add(new Gamete(hs, Tools.seq(g*p2,
                        (g+1)*p2-1), popdata));
                   //each gamete gets ploidy/2 HaploStruct[]
            }
        } else if (popdata.unreducedGametes == 1) {
            //FDR, with normal recombination
            //2 gametes, each formed by combining one normal gamete from the
            //"top" with one from the "bottom"
            //This is equivalent to re-sampling two 2n gametes from
            //the 2*ploidy homologues per chromosome after recombination
            //and chromatid separation,
            //where each gamete must get one copy of the centromere of each
            //parental homologue.
            int[] bottomix = {2, 3}; //index to bottom gametes
            if (rand.nextFloat() < 0.5) bottomix = new int[] {3, 2};
            int[] gam1 = new int[popdata.ploidy];
            int[] gam2 = new int[popdata.ploidy];
            for (int h=0; h<p2; h++) {
                gam1[h] = h; 
                gam1[p2+h] = bottomix[0]*p2+h;
                gam2[h] = p2+h; 
                gam2[p2+h] = bottomix[1]*p2 + h;
            }
            popdata.ploidy = popdata.ploidy * 2;
            gametes.add(new Gamete(hs, gam1, popdata));
            gametes.add(new Gamete(hs, gam2, popdata));
            popdata.ploidy = popdata.ploidy / 2;
        } else if (popdata.unreducedGametes == 2) {
            //SDR, with normal recombination 
            //2 gametes, combining the two normal top or bottom gametes
            //This effectively undoes the second meiotic division
            int[] gam1 = Tools.seq(0, popdata.ploidy-1);
            int[] gam2 = Tools.seq(popdata.ploidy, 2*popdata.ploidy-1);
            popdata.ploidy = popdata.ploidy * 2;
            gametes.add(new Gamete(hs, gam1, popdata));
            gametes.add(new Gamete(hs, gam2, popdata));
            popdata.ploidy = popdata.ploidy / 2;
        }    
        // finally, store the meiotic conficuration in each of the gametes:
        for (Gamete gam: gametes) {
            gam.parconfig = allchromconf;
        }
        //for debugging: print the gametes and recombination positions
        if (false) { //for debugging output: make true
            int[]recCount = new int[popdata.chromCount()];
            int gametehomcount = popdata.ploidy/2;
            if (popdata.unreducedGametes > 0) {
                //with FDR or SDR the gametes have ploidy homologs
                gametehomcount = popdata.ploidy;
            }    
            System.out.println("gam\tchr\thom\tHaploStruct");
            for (int g=0; g<gametes.size(); g++) {
                for (int c=0; c<popdata.chromCount(); c++) {
                    for (int hom=0; hom<gametehomcount; hom++) {
                        System.out.println(g+"\t"+c+"\t"+hom+"\t"+
                        gametes.get(g).getHaploStruct(c, hom).toString());
                        if (gametes.get(g).getHaploStruct(c, hom).isRecombinantHead(popdata.getChrom(c).getTailPos())) {
                            recCount[c]++;
                        }
                    }
                }    
            }
            //list all recombination positions in this meiosis:
            //using a Genotype composed of all 4 (or 2) gametes
            HaploStruct[][] allhs = new HaploStruct[popdata.chromCount()][2*popdata.ploidy];
            for (int c=0; c<popdata.chromCount(); c++) {
                for (int g=0; g<gametes.size(); g++) {
                    for (int hom=0; hom<gametehomcount; hom++) {
                        allhs[c][gametehomcount*g+hom] = gametes.get(g).getHaploStruct(c, hom);
                    }    
                }
            }
            DecimalFormat fix = new DecimalFormat("#0.0000",new DecimalFormatSymbols(Locale.US));
            for (int c=0; c<popdata.chromCount(); c++) {
                double[] recpos = new Genotype(2*popdata.ploidy,allhs,popdata).recPositions(c);
                System.out.print("chrom\t"+c+"\t freqRec=\t"+((2.0*recCount[c])/gametes.size()/popdata.ploidy)+
                        "\trecomb:\t"+recpos.length);
                for (int i=0; i<recpos.length; i++) {
                    System.out.print("\t"+fix.format(recpos[i]));
                }
                System.out.println();
            } 
        }//debug output
        return gametes;
    } //doMeiosis    
    
    /**
     * calcGenomehs makes a 4-dim ArrayList of Integers, with for each chromosome
     * and both sides an array of all ploidy/2 different genome numbers, and
     * for each genome number all the haplostructs
     */
    private void calcGenomehs() {
        genomenr = new int[popdata.chromCount()][][];
        genomehs = new ArrayList<ArrayList<ArrayList<ArrayList<Integer>>>>(); //for each genome nr:
            //which haplostructs in initorder have this genome nr (for both sides) 
        possiblePairCount = new int[popdata.chromCount()][];
        for (int chr=0; chr<popdata.chromCount(); chr++) {
            genomehs.add(new ArrayList<ArrayList<ArrayList<Integer>>>());
            //TetraploidChromosome chrom = (TetraploidChromosome)popdata.getChrom(chr);
            //if (chrom.getPrefPairingProb() > 0.0) {
            genomenr[chr] = new int[2][popdata.ploidy];
            possiblePairCount[chr] = new int[] {0,0}; //the number of possible preferential pairs
                //(for both sides)
            for (int side=0; side<2; side++) {
                genomehs.get(chr).add( new ArrayList<ArrayList<Integer>>());
                for (int g=0; g<popdata.ploidy/2; g++) {
                    genomehs.get(chr).get(side).add(new ArrayList<Integer>());
                }
                for (int h=0; h<popdata.ploidy; h++) {
                    //first get the genome number at the haplostruct side:
                    genomenr[chr][side][h] = (side==0 ?
                            haplostruct[chr][h].getFounderAtHead() :
                            haplostruct[chr][h].getFounderAtTail())
                            % (popdata.ploidy/2);
                    //and list all the haplostruct sides in genomehs:
                    genomehs.get(chr).get(side).get(genomenr[chr][side][h]).add(h);
                }
                //count the number of possible prefpairs at this side:
                for (int g=0; g<popdata.ploidy/2; g++) {
                    possiblePairCount[chr][side] += 
                            genomehs.get(chr).get(side).get(g).size()/2;
                }
            } //for side 
        } //for chr    
    } //calcGenomehs

    /**
     * getPairing returns an int array describing which haplostruct pairs
     * with which other haplostruct at the side (0 or 1, head or tail) of
     * chromosone chr. Pairing is done taking the value of prefPairingProb
     * and the actual configuration of genome numbers into account.
     * @param chr the chromosome number
     * @param side 0 or 1, head or tail
     * @param initorder a randomized order of the haplostructs
     * @return int array of length popdata.ploidy. pairedWith[[i]==j means that 
     * haplostruct[initorder[i]] pairs with haplostruct[initorder[j]]
     */
    private int[] getPairing(int chr, int side) {
        assert(popdata.ploidy > 2);
        int[] pairedWithSide = new int[popdata.ploidy];
        boolean done[] = new boolean[popdata.ploidy];
        TetraploidChromosome tchrom = (TetraploidChromosome)popdata.getChrom(chr);
        if (genomehs == null) calcGenomehs();
        calcGenohs(chr,side);
        int firsth = 0, secondh,  //first and second haplostruct in initorder forming a pair
                gnm, //genome number at chrom [initorder[firsth]], side=side
                hIndex; //index within genohs.get(gnm)
        while (firsth < popdata.ploidy-1) {
            while (firsth < popdata.ploidy-1 && done[firsth]) firsth++;
            if (firsth < popdata.ploidy-1) {
                done[firsth] = true;
                gnm = genomenr[chr][side][initorder[firsth]];
                hIndex = genohs.get(gnm).indexOf(initorder[firsth]);
                genohs.get(gnm).remove(hIndex);
                if ( !genohs.get(gnm).isEmpty() &&
                        rand.nextDouble()<tchrom.getPrefPairingProb() ) {
                    //preferential pairing, get a partner with same genome number:
                    hIndex = rand.nextInt(genohs.get(gnm).size());
                    secondh = backorder[genohs.get(gnm).get(hIndex)];
                } else {
                    //no preferential pairing: get any still available partner
                    secondh = firsth;
                    while (done[secondh]) {
                        secondh = rand.nextInt(popdata.ploidy-firsth-1)+firsth+1;
                    }
                    gnm = genomenr[chr][side][initorder[secondh]];
                    hIndex = genohs.get(gnm).indexOf(initorder[secondh]);
                }
                done[secondh] = true;
                genohs.get(gnm).remove(hIndex);
                pairedWithSide[firsth] = secondh;
                pairedWithSide[secondh] = firsth;
                firsth++;
            }    
        } //while firsth
        //} //preferential pairing    
        return pairedWithSide;
    } //getPairing
    
    /**
     * Pairing is done at head and tail ends of the chromosomes.
     * The pairing at primaryside is always used; the pairing at
     * secondaryside can be overruled if it is not needed or if it
     * conflicts with the specified number of bivalents and quadrivalents.
     * The side where more preferential pairs are possible gets a higher
     * probability to be used as primaryside.
     * Primaryside is determined again for each meiosis.
     * @param chr the chromosome number
     * @return 
     */
    private int getPrimarySide(int chr) {
        double probside0 = possiblePairCount[chr][0] == possiblePairCount[chr][1] ? 
                0.5 :
                (double)(possiblePairCount[chr][0]) /
                    (double)(possiblePairCount[chr][0]+possiblePairCount[chr][1]);
        return rand.nextDouble() < probside0 ? 0 : 1;
    } //getPrimarySide
    
    /* Further variables that need to be global as they are to be accessed 
     * and changed by helper functions of calcChromConfig
     */
    
    /** initorder is a private array of length ploidy. It represents a randomized
     * permutation of the haplostructs of the current chromosome in the current meiosis.
     * E.g. a particular Bivalent might be composed of the haplostructs
     * initorder[firsth] and initorder[secondh].
     */
    private int[] initorder = new int[popdata.ploidy];
    
    /**
     * backorder is the inverse of initorder: initorder[backorder[h]] == h.
     */
    private int[] backorder = new int[popdata.ploidy]; //compiler message incorrect	
    
    private int primaryside, secondaryside; //see getPrimarySide
    private int[][] pairedWith = new int[2][]; //for each haplostruct in initorder:
            //with which other haplostruct is it paired (for both sides)
    private int bivalentcount, quadrivalentcount; //number realized so far
    private ChromConfig chrconf = new ChromConfig(popdata.ploidy); //for the current chromosome
    private boolean[] hsdone; //for each haplostruct in initorder: has it
        //already been used in a Bivalent or Quadrivalent?
    ArrayList<ArrayList<Integer>> genohs;
    
    /**
     * calcGenohs makes a deep copy of one chromosome/side from genomehs
     * and removes all haplostruct entries that have already been used.
     * @param chr the chromosome number
     */
    private void calcGenohs(int chr, int side) {
        genohs = new ArrayList<ArrayList<Integer>>();
        for (int g=0; g<genomehs.get(chr).get(side).size(); g++) {
            genohs.add(new ArrayList<Integer>());
            for (int i=0; i<genomehs.get(chr).get(side).get(g).size(); i++) {
                genohs.get(g).add(genomehs.get(chr).get(side).get(g).get(i));
                //even though we are now referring to the same Integer object as in genomehs
                //it won't matter if that object is later remover from genohs
                //as it will still be present in genomehs
            }
        }
        int gnm, hIndex;
        for (int h=0; h<popdata.ploidy; h++) {
            if (hsdone[h]) {
                gnm = genomenr[chr][secondaryside][initorder[h]];
                hIndex = genohs.get(gnm).indexOf(initorder[h]);
                genohs.get(gnm).remove(hIndex);
            }    
        }
    } //calcGenohs
    
    /**
     * calcChromConfig calculates for one of the (polyploid) chromosomes
     * the ChromConfig
     * @param c the chromosome number
     * @return 
     */
    private void calcChromConfig(int c) {
        chrconf = new ChromConfig(popdata.ploidy);
        chrconf.chromseq = new int[popdata.ploidy];
        if (popdata.ploidy == 2) {
            chrconf.chromseq = tools.shuffleArray(Tools.seq(2));
            chrconf.quad = 0;
        }
        else {
            assert (popdata.ploidy>2);
            TetraploidChromosome chr = (TetraploidChromosome)popdata.getChrom(c);
            initorder = tools.shuffleArray(Tools.seq(popdata.ploidy)); //global but only relevant during
                //this method
            for (int h=0; h<popdata.ploidy; h++) {
                    backorder[initorder[h]] = h;
            }
            hsdone = new boolean[popdata.ploidy];
            for (int side=0; side<2; side++) pairedWith[side] = getPairing(c,side); 
            primaryside = getPrimarySide(c);
            secondaryside = 1 - primaryside;
            bivalentcount = 0;
            quadrivalentcount = 0;

            if (popdata.naturalPairing) {
                //according to the prefpairing of this chromosome
                //we assign first spontaneous Bivalents (i.e. same pairing
                //at primary and secondary side),
                getSpontaneousBivalents(popdata.ploidy/2);
                //Next we search spontaneous Quadrivalents:
                getSpontaneousQuadrivalents((popdata.ploidy-2*bivalentcount)/4);

                //Now we have extracted all spontaneous Bivalents and Quadrivalents.
                //For the remaining haplostruct we ignore the previous pairings at 
                //secondaryside but make new ones ad hoc:
                calcGenohs(c,secondaryside);
                while (2*bivalentcount + 4*quadrivalentcount < popdata.ploidy) {
                    if (2*bivalentcount + 4*quadrivalentcount == popdata.ploidy - 2) {
                        //only 2 chromosomes left
                        getForcedBivalent();
                    } else {
                        //make a new Bivalent or Quadrivalent:
                        getForcedBiOrQuadri(c, true);
                    }    
                } //while (2*bivalentcount + 4*quadrivalentcount < popdata.ploidy) 
            } //if (popdata.naturalPairing)     

            else {
                //no naturalpairing but specified fraction Quadrivalents
                int numQuadrivalents = chr.ranQuadri();
                //now we know the number of quadrivalents to form 
                //with this chromosome in this meiosis
                int numBivalents = (popdata.ploidy-4*numQuadrivalents)/2;
                //we assign first the spontaneous Bivalents and Quadrivalents
                //as far as needed,
                //then add forced Bivalents: use the pairing at primaryside only
                getSpontaneousBivalents(numBivalents);
                getSpontaneousQuadrivalents(numQuadrivalents);
                while (bivalentcount < numBivalents) {
                    getForcedBivalent();
                } 
                if (quadrivalentcount < numQuadrivalents) {
                    calcGenohs(c, secondaryside);
                    while (quadrivalentcount < numQuadrivalents) {
                        getForcedBiOrQuadri(c, false);
                    }
                }
           } //no naturalpairing
            chrconf.quad = quadrivalentcount;
            chr.quadrivalentConfigCount[chrconf.quad]++;
        } //ploidy==2 else
    } //getChrConfig

    private void addBivalent(int firsth, int secondh) {
        //bivalents are at end of chrconf.chromseq
        bivalentcount++;
        chrconf.chromseq[popdata.ploidy-2*bivalentcount] = initorder[firsth];
        chrconf.chromseq[popdata.ploidy-2*bivalentcount+1] = initorder[secondh];
    } //addBivalent
    
    private void addQuadrivalent(int firsth, int secondh, int thirdh, int fourth) {
        /* now the quadrivalent haplostruct are in the order:
            * firsth-secondh-thirdh-fourth
            * with firsth/fourth and secondh/thirdh
            * each paired at primaryside, and
            * firsth/secondh and thirdh/fourth
            * each paired at secondaryside
            * If primaryside==0 (heads) this is ok, but if primaryside==1 (tails)
            * then the haplostruct must be cycled 1 position
            * 
            * Also, the quadrivalents are at the start of chrconf.chromseq
            */
        if (primaryside==0) {
            chrconf.chromseq[(quadrivalentcount*4)] = initorder[firsth];
            chrconf.chromseq[(quadrivalentcount*4)+1] = initorder[secondh];
            chrconf.chromseq[(quadrivalentcount*4)+2] = initorder[thirdh];
            chrconf.chromseq[(quadrivalentcount*4)+3] = initorder[fourth];
        } else {
            chrconf.chromseq[(quadrivalentcount*4)] = initorder[secondh];
            chrconf.chromseq[(quadrivalentcount*4)+1] = initorder[thirdh];
            chrconf.chromseq[(quadrivalentcount*4)+2] = initorder[fourth];
            chrconf.chromseq[(quadrivalentcount*4)+3] = initorder[firsth];
        }  
        quadrivalentcount++;
    } //addQuadrivalent
    
    private void getSpontaneousBivalents(int maxBivalents) {
        int firsth = 0;
        int secondh;
        while (firsth < popdata.ploidy-1 && bivalentcount<maxBivalents) {
            //search next spontaneous bivalent:
            while ( firsth<popdata.ploidy-1 && 
                    ( hsdone[firsth] ||
                        pairedWith[primaryside][firsth] != pairedWith[secondaryside][firsth] ) ) {
                firsth++;
            }
            if (firsth < popdata.ploidy-1) {
                //available firsth and secondh found so bivalent possible
                secondh = pairedWith[primaryside][firsth];
                assert(!hsdone[secondh]);
                addBivalent(firsth, secondh);
                hsdone[firsth] = true;
                hsdone[secondh] = true;
            }
            firsth++;
        } //while firsth, get spontaneous bivalents
    } //getSpontaneousBivalents 
    
    private void getSpontaneousQuadrivalents(int maxQuadrivalents) {
        /* The procedure below assumes that there are no spointaneous bivalents
         * left over. If there is no naturalpairing that may not be the case.
         * Therefore we temporarily mark the spontaneous bivalents in hsdone,
         * then obtain the quadrivalents, and finally put back the
         * temporary bivalents
         */
        ArrayList<Integer> tmpBi = new ArrayList<Integer>();
        int firsth = 0, secondh;
        while (firsth<popdata.ploidy-1) {
            while ( firsth<popdata.ploidy-1 && 
                    ( hsdone[firsth] ||
                        pairedWith[primaryside][firsth] != pairedWith[secondaryside][firsth] ) ) {
                firsth++;
            }
            if (firsth < popdata.ploidy-1) {
                //available firsth and secondh found so bivalent possible
                secondh = pairedWith[primaryside][firsth];
                assert(!hsdone[secondh]);
                hsdone[firsth] = true;
                hsdone[secondh] = true;
                tmpBi.add(firsth);
                tmpBi.add(secondh);
            }
            firsth++;
        }
        // Now we have marked the spontaneous bivalent members in hsdone
        // and continue to find spontaneous quadrivalents:
        firsth = 0;
        int thirdh, fourth;
        while (firsth < popdata.ploidy-3 && quadrivalentcount < maxQuadrivalents) {
            //search next spontaneous Quadrivalent:
            while ( firsth<popdata.ploidy-3 && hsdone[firsth]) firsth++;
            if (firsth < popdata.ploidy-3) {
                //free firsth found so quadrivalent possible
                secondh = pairedWith[secondaryside][firsth];
                assert(!hsdone[secondh]);
                thirdh = pairedWith[primaryside][secondh];
                fourth = pairedWith[primaryside][firsth];
                assert(!hsdone[thirdh]);
                assert(!hsdone[fourth]);
                //These are the four partners in the new potential Quadrivalent.
                //they cannot yet be part of an existing Bivalent or Quadrivalent,
                //else hsdone[firsth] would be true
                if (pairedWith[secondaryside][thirdh] == fourth) {
                    //These 4 form a spontaneous Quadrivalent
                    addQuadrivalent(firsth, secondh, thirdh, fourth);
                    hsdone[firsth] = true;
                    hsdone[secondh] = true;
                    hsdone[thirdh] = true;
                    hsdone[fourth] = true;
                } //spontaneous Quadrivalent is ok
            } //free firsth found so quadrivalent possible    
            firsth++;
        } //while firsth, get spontaneous quadrivalents  
        // Finally we unmark the (extra) spontaneous bivalents in hsdone:
        for (Integer i : tmpBi) {
            assert(hsdone[i]);
            hsdone[i] = false;
        }
    } //getSpontaneousQuadrivalents  

    private void getForcedBivalent() {
        int firsth = 0;
        while (firsth<hsdone.length && hsdone[firsth]) firsth++;
        assert firsth<hsdone.length : "getForcedBivalent";
        int secondh = pairedWith[primaryside][firsth];
        assert(!hsdone[secondh]);
        //next available bivalent pair found
        addBivalent(firsth, secondh);
        hsdone[firsth] = true;
        hsdone[secondh] = true;
    } //getForcedBivalent
    
    /**
     * getForcedBiOrQuadri is called if naturalpairing, after the spontaneous 
     * Bivalents and Quadrivalents have been taken based on pairedWith at
     * both primaryside and secondaryside.
     * The approach is to use the pairedWith pairs at primaryside and to obtain
     * a new, still available partner for firsth at secondaryside taking 
     * the prefPairingProb into account.
     * So:
     * 1. we take the first still available haplostruct -> firsth
     * 2. we take its partner based on pairedWith[promaryside] -> fourth
     * 3. we get a new partner for firsth at secondaryside (not using
     *    pairedWith[secomdaryside] taking preferential pairing into account
     *    -> secondh
     * 4. if secondh==fourth we have a bivalent,  else we take secondh's partner
     *    based on pairedWith[primaryside] -> thirdh and we have a quadrivalent,
     *    imposing pairing between thirdh and fourth at secondaryside
     * @param chr the chromosome number (needed because we have to access
     * genomenr again)
     * @param biAllowed if true either a Bivalent or a Quadrivalent is extracted,
     * if false only a Quadrivalents can be extracted
     */
    private void getForcedBiOrQuadri(int chr, boolean biAllowed) {
        int firsth = 0; 
        int secondh, thirdh, fourth;
        while ( firsth<popdata.ploidy-1 && hsdone[firsth]) firsth++;
        fourth = pairedWith[primaryside][firsth];
        assert(!hsdone[fourth]);
        //we will do preferential pairing at secondary side, if possible, 
        //only with firsth. Therefore if fourth has more options to pair 
        //than firsth we need to switch firsth and fourth:
        int gnm = genomenr[chr][secondaryside][initorder[firsth]];
        int gnm4 = genomenr[chr][secondaryside][initorder[fourth]];
        if (genohs.get(gnm4).size() > genohs.get(gnm).size()) 
         {
            int i=firsth; firsth=fourth; fourth=i;
            i=gnm; gnm=gnm4; gnm4=i;
        }
        int hIndex = genohs.get(gnm).indexOf(initorder[firsth]);
        genohs.get(gnm).remove(hIndex);
        hsdone[firsth] = true;
        int hIndex4 = genohs.get(gnm4).indexOf(initorder[fourth]);
        TetraploidChromosome tchrom = (TetraploidChromosome)popdata.getChrom(chr);
        if ( ( (biAllowed && genohs.get(gnm).size()>0) ||
               (!biAllowed && genohs.get(gnm).size()>1) ) &&
             rand.nextDouble() < tchrom.getPrefPairingProb() ) {
            //preferential pairing:
            do {
                hIndex = rand.nextInt(genohs.get(gnm).size());
                secondh = backorder[genohs.get(gnm).get(hIndex)];
            } while (secondh==fourth && !biAllowed);
            hsdone[fourth] = true;
            genohs.get(gnm4).remove(hIndex4);
            if (secondh == fourth) {
                //new Bivalent created
                addBivalent(firsth, secondh);
            } //prefpairing, new bivalent
            else {
                //No Bivalent: make a Quadrivalent by ignoring the
                //pairing between thirdh and fourth at secondaryside
                assert(!hsdone[secondh]);
                gnm = genomenr[chr][secondaryside][initorder[secondh]];
                hIndex = genohs.get(gnm).indexOf(initorder[secondh]);
                genohs.get(gnm).remove(hIndex);
                hsdone[secondh] = true;
                thirdh = pairedWith[primaryside][secondh];
                assert(!hsdone[thirdh]);
                gnm = genomenr[chr][secondaryside][initorder[thirdh]];
                hIndex = genohs.get(gnm).indexOf(initorder[thirdh]);
                genohs.get(gnm).remove(hIndex);
                hsdone[thirdh] = true;
                addQuadrivalent(firsth, secondh, thirdh, fourth);
            } //prefpairing, new quadrivalent
        } else {
            //no preferential pairing possible: get a random secondh and see 
            //if it is the same as fourth:
            secondh = firsth; 
            while ( hsdone[secondh] ||
                    (secondh==fourth && !biAllowed) ) {
                secondh = rand.nextInt(popdata.ploidy);
            }
            hsdone[fourth] = true;
            genohs.get(gnm4).remove(hIndex4);
            if (secondh==fourth) {
                //new Bivalent found
                addBivalent(firsth, secondh);
                hsdone[secondh] = true;
                //genohs already done
            } //no prefpairing, new bivalent
            else {
                //new Quadrivalent:
                thirdh = pairedWith[primaryside][secondh];
                assert(!hsdone[secondh] && !hsdone[thirdh]);
                hsdone[secondh] = true;
                gnm = genomenr[chr][secondaryside][initorder[secondh]];
                hIndex = genohs.get(gnm).indexOf(initorder[secondh]);
                genohs.get(gnm).remove(hIndex);
                hsdone[thirdh] = true;
                gnm = genomenr[chr][secondaryside][initorder[thirdh]];
                hIndex = genohs.get(gnm).indexOf(initorder[thirdh]);
                genohs.get(gnm).remove(hIndex);
                addQuadrivalent(firsth, secondh, thirdh, fourth);
            } //no prefpairing, new quadrivalent
        } //no prefpairing possible
    } //getForcedBiOrQuadri   
 

}

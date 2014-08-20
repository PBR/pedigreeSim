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
     * prefPairingConfig:
     * used only if popdata.ploidy==4. In that case indicates the situation
     * for each chromosome with respect to preferential pairing.
     * NOTE that if preferential pairing occurs, the even and the odd
     * founder alleles are considered as the two homeologous genomes.
     * int[] with the index referring to the chromosome; contents are:
     * -1 no preferential pairing possible at either end (diploid, or
     *    tetraploid if not one odd and one even pair of founder alleles
     *    at either end)
     * 0  preferential pairing possible at chromosome start only
     * 1  preferential pairing possible at chromosome end only
     * 2  preferential pairing possible at both ends, pairing matches
     *    (at both ends either 01/23 or 02/13, where odd and even may be
     *     the same or reversed in both ends of all chromosomes)
     * 3  preferential pairing possible at both end, pairing does not match
     *    (i.e. at one end the pairing is 01/23, at the other 02/13)
     */
    private int[] prefPairingConfig;

    /**
     * prefPairing:
     * array indicating which chromosome ends may form preferential pairs
     * first index refers to the chromosome
     * second index is 0/1 for chromosome start or end
     * The content is a number from 0 to 3 indicating the preferential pairing:
     * -1 = no pref pairing: 3 or 4 haplostructs have same type (even or odd)
     *  0 = chrom 0/1 vs chrom 2/3
     *  1 = chrom 0/2 vs chrom 1/3
     *  2 = chrom 0/3 vs chrom 1/2
     * (i.e. the content is (partner of chrom 0)-1, or -1 if no pairing)
     * The array always exists
     */
    int[][] prefPairing;

    /**
     * pairs[i-1] list the chromosomes indicated by prefPairing i:
     * the first two numbers are the chromosomes in the first pair,
     * the last two are the chromosomes in the second pair
     */
    static int[][]pairs = new int[][] { {0,1,2,3},{0,2,1,3},{0,3,1,2} };

    /**
     * Class ChromConfig is only relevant for tetraploid individuals.
     * It describes for one chromosome in one meiosis
     * - if the four homologs form a Quadrivalent or two Bivalents
     * - the order of the four homologs:
     *   for Bivalents the first two form the first Bivalent, the second two the second
     *   for a Quadrivalent: the sequence in the cross-type Quadrivalent
     */
    static class ChromConfig {
        boolean quad;
        int[] chromseq;
    }



    /**
     * Constructor to create a founder: its [chromCount][ploidy] HaploStructs
     * each have only one segment, with the given founderAlleles and its parents
     * are set to {null,null}
     * @param name
     * @param popdata
     * @param founderAlleles array of <ploidy> integers
     * @throws Exception 
     */
    public Individual(String name, int[] founderAlleles,
            PopulationData popdata)
            throws Exception {
        super(popdata);
        this.tools = popdata.tools;
        this.rand = popdata.tools.rand;
        this.indivName = name;
        if (founderAlleles.length !=  popdata.ploidy) {
            throw new Exception("Individual constructor: length of founderAlleles incorrect");
        }
        else {
            haplostruct = new HaploStruct[popdata.chromCount()][popdata.ploidy];
            for (int p=0; p<popdata.ploidy; p++) {
                for (int c=0; c<popdata.chromCount(); c++) {
                    try {
                        haplostruct[c][p] = new HaploStruct(popdata.getChrom(c), founderAlleles[p]);
                    } catch (Exception ex) {
                        Logger.getLogger(Individual.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }
            calcPrefPairing();
        }
        this.parents = new Individual[] {null, null};
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
        this(name, haplostruct, popdata);
        setParents(parents);
    }

    /**
     * Constructor specifying the haplostruct but no parents (these are set
     * to {null,null}).
     * This would be the normal constructor to use when reading haplostructs
     * from file, for generating genotypes; the parents are then not
     * relevant.
     * @param name
     * @param haplostruct
     * @param popdata
     * @throws Exception
     */
    public Individual(String name, HaploStruct[][] haplostruct,
        PopulationData popdata)
        throws Exception {
        super(popdata);
        if (name==null || name.isEmpty()) {
            throw new Exception("Individual constructor: name  empty");
        }
        this.indivName = name;
        setHaploStruct(haplostruct);
        calcPrefPairing();
        this.parents = new Individual[] {null,null};
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

    protected void setHaploStruct(HaploStruct[][] haplostruct) throws Exception {
        if (checkHaploStruct(popdata.ploidy, haplostruct)) {
            this.haplostruct = haplostruct;
            calcPrefPairing();
        } else {
            throw new Exception("Individual setHaploStruct: invalid haplostruct");
        }
    }

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

    public boolean isFounder() {
        return parents[0]==null;
    }


    /**
     * mating creates a new genotype from the fusion of two Gametes,
     * the first derived from this parent, the second from otherParent;
     * then thisGamete.fertilization(offspringName, otherGamete) is called.
     * In the offspring the haplostruct from thisGamete are added first
     * and those from otherGamete last; the order of the haplostruct are
     * not changed within or between the gametes.
     * @param otherParent
     * @return
     * @throws Exception
     */
    public HaploStruct[][] mating(Individual otherParent)
            throws Exception {
        if (otherParent==null ||
                !otherParent.getPopdata().equals(popdata)) {
            throw new Exception("Other parent invalid in mating");
        }
        Gamete gamete0 = this.doMeiosis().get(0);
        Gamete gamete1 = otherParent.doMeiosis().get(0);
        return gamete0.fertilization(gamete1);
    }

    /**
     * selfing does the same as mating, with this as otherParent
     * @param offspringName
     * @return
     * @throws Exception
     */
    public HaploStruct[][] selfing(String offspringName) throws Exception {
        Gamete gamete0 = this.doMeiosis().get(0);
        Gamete gamete1 = this.doMeiosis().get(0);
        return gamete0.fertilization(gamete1);
    }

    /**
     * preferential pairing occurs between the telomeres (start and end of
     * chromosome) that both have an odd or both an even founder allele,
     * and only if there are two odd and two even telomeres in the tetraploid
     * individual
     */
    protected void calcPrefPairing() {
        prefPairingConfig = new int[popdata.chromCount()];
        prefPairing = new int[popdata.chromCount()][2];
        for (int chr =0; chr<popdata.chromCount(); chr++) {
            prefPairingConfig[chr] = -1;
            prefPairing[chr] = new int[] {-1,-1};
        }
        if (popdata.ploidy==4) {
            for (int chr=0; chr<popdata.chromCount(); chr++) {
                boolean[][] isOdd = new boolean[4][2]; //preferential pairing between odd and between even founder alleles
                for (int i=0; i<4; i++) {
                    //start, end odd?
                    isOdd[i][0] = haplostruct[chr][i].getFounderAtStart() % 2 == 1;
                    isOdd[i][1] = haplostruct[chr][i].getFounderAtEnd() % 2 == 1;
                }
                for (int side=0; side<2; side++) {
                    //prefPairs[chr][side] = null;
                    int i=1;
                    while (i<4 && isOdd[i][side] != isOdd[0][side]) {
                        i++;
                    }
                    if (i<4) {
                        //else no pairing for haplostruct 0
                        //now 0 and i are the first pair
                        int j= i==1 ? 2 : 1;
                        if (isOdd[j][side] != isOdd[0][side]) {//else at leat three odd or even
                            int k = (i==2 || j==2) ? 3 : 2;
                            if (isOdd[k][side] == isOdd[j][side]) { //else no two different pairs
                                //0 and i are the first pair, j and k the second;
                                //the pairs and the homologues within the pairs
                                //are ordered, there are only 3 possibilities:
                                //{0,1,2,3}, {0,2,1,3}, {0,3,1,2}
                                prefPairing[chr][side] = i-1; //chrom 0/i vs chrom j/k

                            }
                        }
                    }
                } //for side
                //determine prefPairingConfig for chr:
                if (prefPairing[chr][0]==-1) {
                    if (prefPairing[chr][1]==-1) prefPairingConfig[chr] = -1;
                    else prefPairingConfig[chr] = 1; //only at side 1
                }
                else {
                    if (prefPairing[chr][1]==-1) prefPairingConfig[chr] = 0; //only at side 0
                    else { //both ends have pref. pairs; matching or not?
                        if (prefPairing[chr][0]==prefPairing[chr][1])
                            prefPairingConfig[chr] = 2;
                        else prefPairingConfig[chr] = 3;
                    }
                }
            } //for chr
        } //ploidy==4
    } //calcPrefPairing

    /**
     * doNaturalPairing performs preferential or random pairing of both ends
     * of chromosome chrom, and based on the result determines if they form
     * two Bivalents or one Quadrivalent and
     * how the Bivalents or Quadrivalent are constructed
     * @param chrom
     * @return a ChromConfig:
     *   if two bivalents, the first two are the first Bivalent and the
     *   second two the second Bivalent;
     *   if one Quadrivalent, the order is such that the correct pairing
     *   of starts and ends occurs where result[0] goes to Quadrivalent slots 0&1,
     *   return[1] to slots2&3 etc
     */
    ChromConfig doNaturalPairing(int chrom) {
        ChromConfig result = new ChromConfig();
        int[] sidePairing = new int[2]; //any of the 3 possible pairings, as in prefPairing, except -1
        for (int side=0; side<2; side++) {
            if ((prefPairing[chrom][side]==-1) ||
                (rand.nextDouble() >
                ((TetraploidChromosome)(popdata.getChrom(chrom))).getPrefPairingProb())) {
                //random pairing at this end
                sidePairing[side] = rand.nextInt(3);
            }
            else {
                //use pairing:
                sidePairing[side] = prefPairing[chrom][side];
            }
            /*now we need to see if the pairs at both ends match or not;
             *at both sides the pairing can be (ordered!)
             * 01/23, 02/13, 03/12
             * if the pairs at start and end are identical we have two bivalents
             * else we have a quadrivalent as follows:
             * start end   slot01 slot23 slot45 slot67
             * 01/23 02/13   0       2     3       1
             * 01/23 03/12   0       3     2       1
             * 02/13 01/23   0       1     3       2
             * 02/13 03/12   0       3     1       2
             * 03/12 01/23   0       1     2       3
             * 03/12 02/13   0       2     1       3
             */
        }
        //do the pairings at the two ends match?
        result.quad = sidePairing[0]!=sidePairing[1];
        if (result.quad) {
            //put chrom in order of slots
            result.chromseq = twoEndOrdersToQuadrivalentOrder(sidePairing);
        }
        else { //two bivalents
            /* the pairings at both ends match so we can use either to define
             * the order; use the first = order[0]:
             */
            result.chromseq = randomizeBivalents(pairs[sidePairing[0]].clone());
        }
        return result;
    } //doNaturalPairing

    /**
     * twoEndOrdersToQuadrivalentOrder
     * This method takes the pairwise ordering of the chromosome homolog
     * starts (order[0] and ends (order[1]) and from it derives the
     * sequence of the chromosomes in a Quadrivalent.
     * Both order[side] are ordered and have one of only three
     * possible sequences: {0,1,2,3}, {0,2,1,3}, {0,3,1,2} (the sequence
     * indicates the first and the second pairs of haplostructs at that side)
     * if the pairs at start and end are identical we have two bivalents
     * and that is not allowed here (see assertion)
     * else we have a quadrivalent as follows:
     * start end   slot01 slot23 slot45 slot67
     * 01/23 02/13   0      2      3      1
     * 01/23 03/12   0      3      2      1
     * 02/13 01/23   0      1      3      2
     * 02/13 03/12   0      3      1      2
     * 03/12 01/23   0      1      2      3
     * 03/12 02/13   0      2      1      3
     * @param order the pairing at both sides of the chromosomes
     * @return the order of the chromosomes in the Quadrivalent (randomized)
     */
    int[] twoEndOrdersToQuadrivalentOrder(int[] sidePairing) {
        assert(sidePairing.length==2);
        assert(sidePairing[0]!=sidePairing[1]); //else the two pairings match and two bivalents would result
        //put chrom in order of slots
        int[] seq = new int[4];
        seq[0] = 0;
        /*seq[1] = order[1][1];
        seq[3] = order[0][1];
        seq[2] = 6 - seq[1] - seq[3];*/
        seq[1] = sidePairing[1]+1; //the chrom paired with chrom 0 at the end
        seq[3] = sidePairing[0]+1; //the chrom paired with chrom 0 at the start
        seq[2] = 6 - seq[1] - seq[3]; //the remaining chrom
        /* now we still must randomize while keeping the paired ends together:
         * we have abcd, the other possibilities with same pairing are:
         * badc : horizontal flip
         * dcba : vertical flip
         * cdab : both flips = rotation 180 deg
         * note that rotation over 90 degrees would exchange starts and ends
         */
        int[] result = new int[4];
        switch (rand.nextInt(4)) {
            case 0 : result = seq; break;
            case 1 : result = new int[] {seq[1],seq[0],seq[3],seq[2]}; break;
            case 2 : result = new int[] {seq[3],seq[2],seq[1],seq[0]}; break;
            default : result = new int[] {seq[2],seq[3],seq[0],seq[1]};
        }
        return result;
    } //twoEndOrdersToQuadrivalentOrder

    /**
     * makeBivalents produces two Bivalents based on preferential or
     * random pairing of chromosome chrom
     * @param chrom
     * @return
     */
    int[] makeBivalents(int chrom) {
        int[] result;
        if (prefPairingConfig[chrom]==-1 || //no pairs at either end
            prefPairingConfig[chrom]==3 ||  //non-matching pairs at the two ends
            (rand.nextDouble() >
                ((TetraploidChromosome)(popdata.getChrom(chrom))).getPrefPairingProb())) {
            //no preferential pairing:
            result = pairs[rand.nextInt(3)].clone();
        }
        else {
            //preferential pairing:
            if (prefPairingConfig[chrom] >= 1) {
                //pairing based on chromosome end
                result = pairs[prefPairing[chrom][1]].clone();
            }
            else {
                //pairing based on chromosome start
                result = pairs[prefPairing[chrom][0]].clone();
            }
        }
        //now randomizing between and within pairs is still needed
        //as in doNaturalPairing:
        return randomizeBivalents(result);
    } //makeBivalents
    
    int[] randomizeBivalents(int[] pairing) {
        /*
         * pairing has the pairing of the four chromosomes
         * into two Bivalents, after applying preferential or
         * non-preferential pairing.
         * now we must still randomize while keeping the pairing intact:
         * the first and second pair can be switched, and within each pair
         * the chromosomes may be switched:
         */
        int tmp;
        if (rand.nextBoolean()) { //switch the two pairs:
            tmp=pairing[0]; pairing[0]=pairing[2]; pairing[2]=tmp;
            tmp=pairing[1]; pairing[1]=pairing[3]; pairing[3]=tmp;
        }
        if (rand.nextBoolean()) { //switch within first pair
            tmp=pairing[0]; pairing[0]=pairing[1]; pairing[1]=tmp;
        }
        if (rand.nextBoolean()) { //switch within second pair
            tmp=pairing[2]; pairing[2]=pairing[3]; pairing[3]=tmp;
        }
        return pairing;
    } //randomizeBivalents

    /**
     * makeQuadrivalent produces one Quadrivalent based on preferential or
     * random pairing of chromosome chrom
     * @param chrom
     * @return
     */
    int[] makeQuadrivalent(int chrom) {
        /* if there are no pairs at either end, first at the start random
         * pairing is applied and then at the end such that the pairing at the
         * end is different from the start
         * if there is pairing at only one end, preferential pairing
         * is applied there first and at the other end a random pairing different
         * from the first pairing
         * if there are matching pairs at both ends, first preferential
         * pairing is done at both ends. If this results in two matching
         * pairings there are three possibilities:
         * - at no one end the actual pairing matches the pairs; then one end
         * is selected randomly and there preferential pairing is applied
         * until non-matching pairing is achieved
         * - at one end the actual pairing matches the pairs; then at the
         * other end random pairing different from the first pairing
         * is done
         * - at both ends the actual pairing matches the pairs; then at one
         * end (randomly selected) random pairing different from the first
         * pairing is done
         *
         */
        //int[][] order = new int[2][];
        int[] sidePairing = new int[2];
        int y;
        int side;
        switch (prefPairingConfig[chrom]) {
            case -1: { //no pref. pairing possible at either end
                sidePairing[0] = rand.nextInt(3); //any of the 3 possible pairings, as in prefPairing, except -1
                y = rand.nextInt(2); //any of 2 pairings left over :
                sidePairing[1] = (sidePairing[0]+y+1) % 3;
                break;
            }
            case 0:
            case 1:
            case 2:  { //only pref.pairing possible at one end
                side = prefPairingConfig[chrom];
                if (side==2) side = rand.nextInt(2); //the two side pairings are conflicting, select one
                if (rand.nextDouble() >
                    ((TetraploidChromosome)(popdata.getChrom(chrom))).getPrefPairingProb()) {
                    //random pairing at this end
                    sidePairing[side] = rand.nextInt(3);
                }
                else {
                    //use pref. pairing:
                    sidePairing[side] = prefPairing[chrom][side];
                }
                //for the other end use random pairing:
                y = rand.nextInt(2); //any of 2 pairings left over :
                sidePairing[1-side] = (sidePairing[side]+y+1) % 3;
                break;
            }
            case 3 : { //pref.pairing possible at both ends, not conflicting
                for (side=0; side<2; side++) {
                    if (rand.nextDouble() >
                        ((TetraploidChromosome)(popdata.getChrom(chrom))).getPrefPairingProb()) {
                        //random pairing at this end
                        sidePairing[side] = rand.nextInt(3);
                    }
                    else {
                        //use pref. pairing:
                        sidePairing[side] = prefPairing[chrom][side];
                    }
                }
                //do the actual pairings at the two ends match?
                if (sidePairing[0]==sidePairing[1]) {
                    //we pick one of the sides and do a random pairing there
                    side = rand.nextInt(2); //side 0 or 1
                    y = rand.nextInt(2); ////any of 2 pairings left over :
                    sidePairing[side] = (sidePairing[1-side]+y+1) % 3;
                }
                break;
            }
            default:
                assert (false) : prefPairingConfig[chrom] ;
        } //switch (prefPairingConfig[chrom])
        int[] result = twoEndOrdersToQuadrivalentOrder(sidePairing);
        return result;
    } //makeQuadrivalent


    /**
     * doMeiosis produces four gametes.
     * If ploidy==2, each gamete contains one HaploStruct per chromosome;
     * if ploidy==4, each gamete contains two HaploStruct per chromosome
     * @return an ArrayList of size 4: 4 Gametes.
     *   each element is a HaploStruct[][], with
     *   first index is chromosome number, and
     *   second index is number of HaploStruct (0 for haploid gamete,
     *   0 and 1 for diploid gamete)
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
        HaploStruct[][] hs = new HaploStruct[popdata.chromCount()][];
        ArrayList<Gamete> gametes = new ArrayList<Gamete>();
        HaploStruct[] origChrom;
            //will contain the chromosomes to go into the Bivalent(s) or Quadrivalent
        if (popdata.ploidy==2) {
            origChrom = new HaploStruct[2];
            for (int c=0; c<popdata.chromCount(); c++) {
                //randomize order of the two haplostruct for the Bivalent:
                int first = rand.nextInt(2);
                origChrom[0] = haplostruct[c][first];
                origChrom[1] = haplostruct[c][1-first];
                //do the meiosis, returning 4 haplostruct (1 per gamete):
                hs[c] = new Bivalent(origChrom).doMeiosis();
            }
            for (int g=0; g<4; g++) {
                gametes.add(new Gamete(hs, new int[]{g}, popdata));
                   //each gamete gets one of the four HaploStruct[] per c(hromosome)
            }
        } else {
            assert (popdata.ploidy==4);
            ChromConfig chromconfig;
            for (int c=0; c<popdata.chromCount(); c++) {
                TetraploidChromosome chr = (TetraploidChromosome) popdata.getChrom(c);
                if (popdata.naturalPairing) {
                    chromconfig = doNaturalPairing(c);
                }
                else {
                    chromconfig = new ChromConfig();
                    // no natural pairing; decide between Quadrivalent or two Bivalents
                    chromconfig.quad = chr.getFracQuadrivalents()==0.0 ? false :
                            chr.getFracQuadrivalents()==1.0 ? true :
                            rand.nextDouble()<chr.getFracQuadrivalents();
                    if (chromconfig.quad) {
                        chromconfig.chromseq = makeQuadrivalent(c);
                    }
                    else {
                        chromconfig.chromseq = makeBivalents(c);
                    }
                }
                if (chromconfig.quad) {
                    origChrom = new HaploStruct[4];
                    for (int i=0; i<4; i++) {
                        origChrom[i] = haplostruct[c][chromconfig.chromseq[i]];
                    }
                    hs[c] = new Quadrivalent(origChrom).doMeiosis();
                }
                else {
                    hs[c] = new HaploStruct[2*popdata.ploidy];
                    //make two Bivalents
                    HaploStruct[] chrgam; 
                    for (int p=0; p<popdata.ploidy/2; p++) {
                        chrgam = new Bivalent(
                                haplostruct[c][chromconfig.chromseq[2*p]],
                                haplostruct[c][chromconfig.chromseq[2*p+1]]
                                ).doMeiosis();
                        /* chrgam has the four haploid gametes for one of the
                         * two bivalents for this chromosome.
                         * Now place them into the hs array:
                         */
                        for (int g=0; g<4; g++) {
                            hs[c][2*g+p] = chrgam[g];
                        }
                    } //for p
                } //two Bivalents
                //next we shuffle the two chromosomes per diploid gamete:
                HaploStruct tmp;
                for (int g=0; g<4; g++) {
                    if (rand.nextBoolean()) {
                        tmp = hs[c][2*g];
                        hs[c][2*g] = hs[c][2*g+1];
                        hs[c][2*g+1] = tmp;
                    }
                }
            } //for c
            //now we have the four gametes represented in hs as 8 HaploStructs per
            //chromosome; these are now converted into Gametes:
            for (int g=0; g<4; g++) {
                gametes.add(new Gamete(hs, new int[]{2*g, 2*g+1}, popdata));
                   //each gamete gets two of the eight HaploStruct[]
            }
        } //ploidy==4
        //for debugging: print the four gametes and chiasma positions
        if (false) {
            int[]recCount = new int[popdata.chromCount()];
            System.out.println("gam\tchr\thom\tHaploStruct");
            for (int g=0; g<gametes.size(); g++) {
                for (int c=0; c<popdata.chromCount(); c++) {
                    for (int hom=0; hom<popdata.ploidy/2; hom++) {
                        System.out.println(g+"\t"+c+"\t"+hom+"\t"+
                        gametes.get(g).getHaploStruct(c, hom).toString());
                        if (gametes.get(g).getHaploStruct(c, hom).isRecombinantStart(popdata.getChrom(c).getEndPos())) {
                            recCount[c]++;
                        }
                    }
                }    
            }
            //list all chiasma positions in this meiosis:
            //using a Genotype composed of all 4 gametes
            HaploStruct[][] allhs = new HaploStruct[popdata.chromCount()][2*popdata.ploidy];
            for (int c=0; c<popdata.chromCount(); c++) {
                for (int g=0; g<gametes.size(); g++) {
                    for (int hom=0; hom<popdata.ploidy/2; hom++) {
                        allhs[c][popdata.ploidy*g/2+hom] = gametes.get(g).getHaploStruct(c, hom);
                    }    
                }
            }
            DecimalFormat fix = new DecimalFormat("#0.0000",new DecimalFormatSymbols(Locale.US));
            for (int c=0; c<popdata.chromCount(); c++) {
                double[] recpos = new Genotype(2*popdata.ploidy,allhs,popdata).recPositions(c);
                System.out.print("chrom\t"+c+"\t freqRec=\t"+((2.0*recCount[c])/gametes.size()/popdata.ploidy)+
                        "\tchiasmata:\t"+recpos.length);
                for (int i=0; i<recpos.length; i++) {
                    System.out.print("\t"+fix.format(recpos[i]));
                }
                System.out.println();
            } 
        }
        return gametes;
    } //doMeiosis

}

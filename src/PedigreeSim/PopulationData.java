/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.util.ArrayList;

/**
 * PopulationData holds most of the global variables in the current
 * population, including the individuals, chromosomes, ploidy and
 * a variety of options affecting the meiosis.
 * Most of these settings are intentionally not private but have
 * default visibility, meaning that they can be accesses from
 * Multivalents, Chromosomes and Genotypes (and their descendant classes
 * in this package) without overheads for setters and getters.
 * @author Roeland Voorrips
 */
public class PopulationData {
    
    int ploidy; //2 or 4
    long randomSeed;
    String missing;

    /**
     * chiasmaInterference:
     * if false, no chiasma interference, leading to Haldane map function
     * in Bivalents. If true, chiasma interference is implemented such that
     * the Kosambi map function applies in Bivalents.
     */
    boolean chiasmaInterference; 

    /**
     * naturalPairing:
     * if true (default),in a tetraploid, decide whether two Bivalents or one
     * Quadrivalent is formed based on the pairing of the telomeres.
     * if false, the ratio of Bivalents to Quadrivalents is determined by
     * each TetraploidChromosome.fradQuadrivalents
     */
    boolean naturalPairing;

    /**
     * parallelMultivalents:
     * By default we use an "arm-based" pairing model resulting in a cross-type
     * quadrivalent with four arms where at each point 2 chromosomes are paired.
     * This parameter specifies the probability of another multivalent 
     * configuration, where all chromosomes are parallel overt their full length
     * and each chiasma can occur between any pair of chromosomes;
     * probably less realistic than cross-type multivalents
     */
    double parallelMultivalents; //default 0.0

    /**
     * pairedCentromeres:
     * applies only to Quadrivalents
     * In a cross-type quadrivalent, two models are suggested in literature 
     * about the way that the four centromeres (and attached, possibly 
     * recombined chromatids) are distributed over the two poles in the first 
     * meiotic division. One model is that the four centromeres can be divided 
     * into two pairs in any of the three possible ways (ab/cd, ac/bd, ad/bc) 
     * with equal probability; the second model assumes that the centromeres 
     * are paired according to the quadrivalent arm in which they end up, 
     * and that from each of the two pairs one centromere moves to either pole. 
     * The second models leads to higher expected amount of double reduction 
     * than the first: in the second model ½  half of the divisions are in the 
     * “adjacent” configuration (the only configuration in which double reduction 
     * can occur), while in the first model only ⅓ of the divisions are 
     * “adjacent”. This parameter specifies the probability of the second model 
     * (range 0.0-1.0); the default is 0.0, i.e. all cross-type quadrivalents 
     * segregate according to the first model.
     */
    double pairedCentromeres; //default: 0.0

    /**
     * bivalentsBidirectional:
     * if false, chiasmata in bivalents are generated from the head
     * to the tail of the chromosome; this is the default situation.
     * if true, they are generated alternating from both ends. This is only
     * an option in the TEST configuration, to check if this approach
     * (which is used for quadrivalents) yields the same distribution
     * of chiasmata positions.
     */
    boolean bivalentsBidirectional; //default: false

    /**
     * allowNoChiasmata:
     * if false, any multivalent meiosis in which one or more chromosomes
     * are not involved in at least one chiasma is discarded and performed
     * again, until all chromosomes are involved.
     * This increases the overall rate of recombination: for bivalents
     * over the whole length of the chromosome the amount of recombination
     * will always be 50% (i.e. the map length becomes infinite).
     * An approach that ensures that the average number of
     * chiasmata in this case will still be 2*chromosome.length (as when no
     * bivalents are rejected) where no chiasma interference takes place
     * is to have m = 0.5/(1-exp(-L/m)), where m is the average chiasma distance
     * and L the chromosome length in Morgan.
     * This is solved numerically
     * in Tools.calcChiasmaDist. But as remarked above, the effective total
     * chromosome map length will become infinite, so over longer distances
     * the marker distances will be estimated larger than actual:
     * for L = 1.0 the  % overestimation is as follows:
     * cm    1   2   5   10  20  30  40   50   60   70   80   90
     * %over 2.3 2.4 2.4 2.5 5.0 7.8 11.9 17.0 23.8 32.9 46.8 73.6
     */
    boolean allowNoChiasmata; //default: true

    /**
     * quadriEachArm:
     * applies only to Quadrivalents
     * if true, chiasmata are generated in turn from the end of each arm
     * (each cycle we get a random order of the 4 arms and generate the
     * chiasmata in that order; we repat until no further chiasmata fit)
     * if false, each new chiasma is generated from a random arm, so several
     * consecutive chiasmata may be generated in the same arm. This is expected
     * to result in a wider fluctuation of the chromosome exchange point
     * (the center of the quadrivalent cross)
     * false would seem the more natural situation and is the default
     */
    boolean quadriEachArm; //default: false
    
    /**
     * quadriMethod:
     * explained in Quadrivalent - doCrossingOver.
     * Briefly:
     * 1 - crossing-over continues until in all arms a conflicting chiasma
     *     has been produced; exchange interval is between the innermost
     *     chiasmata.
     * 2 - crossing-over continues until in all arms a conflicting chiasma
     *     has been produced. Difference with 1 is that a failed chiasma
     *     does narrow the remaining exchange interval, so that in the 
     *     remaining arms the room for valid chiasmata decreases.
     * 3 - crossing-over stops when the first conflicting chiasma is
     *     produced. The exchange interval is between the innermost
     *     chiasmata.
     */
    int quadriMethod;
    
    int founderAlleleCount; //is ploidy * number of founders
    private ArrayList<Chromosome> chromosome;
    private ArrayList<Individual> individual;
    public Tools tools;

    boolean testMode; //true: test mode, false: normal mode
    //variables only used in test mode:
    int testIter, //onumber of iterations (of testMeioseCount meioses)
        testMeioseCount, //number of meioses per iteration
        testPrintGametes; //0=no, 1=only one (first), else all gametes   
    boolean 
        testPrintEachIter, //print cumulative statistics for each iteration
        testPrintMapdistances, //print map distances in addition to recombination freqs
        testPrintPubTables; //only in combination with specific chrom and map files   
    
    /**
     * The full constructor
     * @param ploidy 2 or 4 (will possibly be extended)
     * @param chiasmaInterference false or true for Haldane or Kosambi map
     * function (in Bivalents; in quadrivalents the recombination fractions
     * and distances are larger)
     * @param parallelMultivalents quadrivalents (and higher) not arm-based
     * but all chromosomes parallel (not realistic)
     * @param allowNoChiasmata if false multivalents in which chromosomes
     * not involved in at least one chiasma ("unpaired chromosomes") are
     * rejected; leads to serious departures from set map distances
     * @param randomSeed allows to re-run same simulation
     * @param missing string to represent missing parents in pedigree
     * and non-existing recombination positions in HaploStruct files
     * @throws Exception
     */
    public PopulationData(int ploidy, boolean chiasmaInterference,
            double parallelMultivalents, double pairedCentromeres,
            boolean naturalPairing,
            boolean allowNoChiasmata,
            boolean bivalentsBidirectional,
            long randomSeed, String missing)
            throws Exception {
        if //(ploidy!=2 && ploidy!=4) {
           (ploidy <= 0 || ploidy % 2 != 0) {     
            throw new Exception("Error in PopulationData: ploidy must be even");
        }
        this.ploidy = ploidy;
        this.chiasmaInterference = chiasmaInterference;
        this.parallelMultivalents = parallelMultivalents;
        this.pairedCentromeres = pairedCentromeres;
        this.bivalentsBidirectional = bivalentsBidirectional;
        this.allowNoChiasmata = allowNoChiasmata;
        this.naturalPairing = naturalPairing;
        this.missing = missing;
        this.randomSeed = randomSeed;
        this.quadriEachArm = false; //default, change only for test purposes
        this.quadriMethod = 2; //2 is default, most realistic
        //set variables for test mode to defaults:
        this.testMode = false;
        this.testPrintGametes = 0;
        this.testPrintEachIter = false;
        this.testPrintMapdistances = true;
        this.testPrintPubTables = false;
        tools = new Tools(randomSeed);
        chromosome = new ArrayList<Chromosome>();
        individual = new ArrayList<Individual>();
    }

    public ArrayList<Chromosome> getChromosome() {
        return chromosome;
    }

    public Chromosome getChrom(int chr) {
        return chromosome.get(chr);
    }

    public int chromCount() {
        return chromosome.size();
    }

    public ArrayList<Individual> getIndividual() {
        return individual;
    }

    public Individual getIndiv(int ind) {
        return individual.get(ind);
    }

    public Individual getIndiv(String name) {
        int i = 0;
        while (i<individual.size() &&
                !individual.get(i).getIndivName().equals(name)) {
            i++;
        }
        if (i>=individual.size()) return null;
        else return individual.get(i);
    }

    public int indivCount() {
        return individual.size();
    }

    public void addChromosome(Chromosome chrom) throws Exception {
        if (chrom==null) {
            throw new Exception("Error in addChromosome: chrom==null");
        }
        if ((ploidy==2 && chrom.getClass().equals(TetraploidChromosome.class)) ||
            (ploidy>2 && !chrom.getClass().equals(TetraploidChromosome.class))) {
            throw new Exception("Error in addChromosome: class of chrom does not match ploidy");
        }
        if (!this.equals(chrom.getPopdata())) {
            throw new Exception("Error in addChromosome: popdata invalid");
        }
        chromosome.add(chrom);
        //note that the same chromosome may be added multiple times, but
        //each instance will be non-homologous to the others
    }

    public void addIndividual(Individual indiv) throws Exception {
        if (indiv==null) {
            throw new Exception("Error in addIndividual: indiv==null");
        }
        if (!indiv.getPopdata().equals(this)) {
            throw new Exception("Error in addIndividual: popdata invalid");
        }
        individual.add(indiv);
    }

    private boolean indPresent(String ind, String[][] indarr, int i) {
        //true if ind occurs in indarr before [i]
        //(by setting i to IndivCount the whole Ind array is checked)
        if (ind==null || ind.equals(missing)) return true;
        int j=i-1;
        boolean found = false;
        while (j>=0 && !found) {
            found = ind.equals(indarr[j][0]);
            j--;
        }
        return found;
    }

    /**
     * sortPedigree takes as input ped, an ArrayList with elements consisting
     * of 3 strings: an individual's name and the names of its 2 parents.
     * The array is sorted such that an individual always appears after its
     * parents, and is then converted to the ArrayList individuals,
     * with the founders receiving founder HaploStructs and all progeny
     * receiving null HaploStructs
     * @param ped ArrayList with elements consisting
     * of 3 strings: an individual's name and the names of its 2 parents
     * @return Error message, "" if no error
     * @throws Exception
     */
    public ArrayList<String[]> sortPedigree(ArrayList<String[]> ped) throws Exception {
        //check if each element consists of 3 Strings:
        int iCount = ped.size();
        int j = 0;
        while (j < iCount && ped.get(j).length==3) {
            j++;
        }
        if (j < iCount) {
            throw new Exception(Integer.toString(j)+"-th line of pedigree does not contain 3 items");
        }
        //check if any indiv is it's own parent:
        j = 0;
        while (j < iCount
                && !ped.get(j)[1].equals(ped.get(j)[0])
                && !ped.get(j)[2].equals(ped.get(j)[0]) ) {
            j++;
        }
        if (j < iCount) {
            throw new Exception("Individual '" + ped.get(j)[0] + " is it's own parent!");
        }
        //check if any indiv has missing name:
        j = 0;
        while (j < iCount
                && !ped.get(j)[0].equals(missing) ) {
            j++;
        }
        if (j < iCount) {
           throw new Exception("Individual '" + j + " has a missing name, which is not allowed");
        }
        //check if any indiv has one but not two parents:
        j = 0;
        while (j < iCount
                && !(ped.get(j)[1].equals(missing) ^
                     ped.get(j)[2].equals(missing) ) ) {
            j++;
        }
        if (j < iCount) {
            throw new Exception("Individual " + ped.get(j)[0] +
                    " has one known and one unknown parent, which is not allowed");
        }
        //check if any individual occurs more than once:
        j = 1; //not 0 !
        while (j < iCount) {
            int k = j-1;
            while (k>=0) {
                if (ped.get(j)[0].equals(ped.get(k)[0])) {
                    throw new Exception("Individual "  + ped.get(j)[0] +
                            " occurs more than once in pedigree");
                }
                k--;
            }
            j++;
        }

        //next: prepare for the sorting:
        String[][] ind = new String[2 * iCount][]; //room at end needed for sorting
        for (int i = 0; i < iCount; i++) {
            ind[i] = ped.get(i);
        }
        /* Sort the individual array in indivCount-1 passes. In each pass i,
         * an individual[i] is placed of which both parents (if not unknown)
         * are already placed before.
         * This individual is searched from position i down. If it is found at
         * position j, all the Indivs at position i to j-1 are moved to the end
         * of the array, and the gap is filled by moving all Indivs j-i
         * positions up.
         * In this way, the order of e.g. sibs placed as one consecutive block
         * is conserved.
         */
        for (int pass=0; pass<iCount; pass++) {
            //for sorting, pass<ICount-1 would be sufficient, but in that case
            //the last Indiv might have a non-existing parent which would not be noticed
            int i = pass;
            //find first Indiv i that can be placed at position pass:
            while ( i<iCount &&
                  ( !indPresent(ind[i][1],ind,i) ||
                    !indPresent(ind[i][2],ind,i) ) ) {
                i++;
            }
            if (i>=iCount) {
                //no Indiv remains with both parents already placed;
                //report the first where one parent lacks:
                i = pass;
                while ( i < iCount &&
                        indPresent(ind[i][1], ind, iCount) &&
                        indPresent(ind[i][2], ind, iCount)) {
                    i++; 
                }
                if (i < iCount) {
                    if (!indPresent(ind[i][1],ind,iCount)) {
                        throw new Exception(ind[i][1] +
                                " is listed as parent of " + ind[i][0] +
                                " but it does not occur in the population");
                    } else {
                        throw new Exception(ind[i][2] +
                                " is listed as parent of " + ind[i][0] +
                                " but it does not occur in the population");
                    }
                } else {
                    throw new Exception("The pedigree could not be sorted, "+
                            "probably due to circular references");
                }
            }
            //if we arrive here without Exception all is well
            if (i>pass) { //else i==pass, ind[i] already correcly placed
                //move all indivs from pass to i-1 to end:
                for (int k=pass; k<i; k++) {
                    ind[iCount+k-pass] = ind[k];
                }
                //move all Indivs from i to new end up to close gap:
                for (int k=0; k<iCount-pass; k++) {
                    ind[pass+k]=ind[i+k];
                }
            }
        } //for pass
        //finished sorting, clear and fill ArrayList individual:
        ArrayList<String[]> result = new ArrayList<String[]>();
        for (int i=0; i<iCount;i++) {
            result.add(ind[i]);
        }
        return result;
    } ////SortPedigree  Compiler warning about CommentChar never used incorrect!

    public String createPopulation(ArrayList<String[]> ped) {
        individual = new ArrayList<Individual>();
        founderAlleleCount=0;
        int iCount = ped.size();
        for (int i=0; i<iCount; i++) {
            //System.out.println(i+": "+ped.get(i)[0]+" "+ped.get(i)[1]+" "+ped.get(i)[2]);
            if (ped.get(i)[1].equals(missing)) {
                //add founder with founderalleles
                int[] fAll = new int[ploidy];
                for (int p=0; p<ploidy; p++) {
                    fAll[p] = founderAlleleCount++;
                }
                try { individual.add(new Individual(ped.get(i)[0], fAll, this));
                } catch (Exception ex) { return ex.getMessage(); }
            }
            else {
                //add non-founder with parents but no HaploStruct
                try {
                    Individual parent1;
                    Individual parent2;
                    try {
                        parent1 = getIndiv(ped.get(i)[1]);
                        parent2 = getIndiv(ped.get(i)[2]);
                    }
                    catch (Exception ex) {
                        throw new Exception("Parent(s) of individual "+ped.get(i)[0]+" not found");
                    }
                    individual.add(new Individual(ped.get(i)[0],
                        new Individual[] {parent1, parent2},
                        this));
                } catch (Exception ex) { return ex.getMessage(); }
            }
        }
        return "";
    } //createPopulation

    public String[] popdataSettings() {
        String[] result = new String[] {
            "ploidy=\t"+ploidy,
            "chiasmaInterference=\t"+chiasmaInterference,
            "testMeioseCount=\t"+testMeioseCount,
            "testIter=\t"+testIter,
            "founderAlleleCount=\t"+founderAlleleCount,
            "missing=\t"+missing,
            "parallelMultivalents=\t"+parallelMultivalents,
            "pairedCentromeres=\t"+pairedCentromeres,
            "quadriMethod=\t"+quadriMethod,
            "allowNoChiasmata=\t"+allowNoChiasmata,
            "naturalPairing=\t"+naturalPairing,
            "bivalentsBidirectional=\t"+bivalentsBidirectional,
            "randomSeed=\t"+randomSeed };
        return result;
    } //popdataSettings

    public int getMaxAlleleCount() {
        int result = 0;
        for (Chromosome chrom : getChromosome()) {
            if (chrom.getLocus() != null) for (Locus loc: chrom.getLocus()) {
                if (loc.getAlleleCount() > result) {
                    result = loc.getAlleleCount();
                }
            }
        }
        return result;
    }

}

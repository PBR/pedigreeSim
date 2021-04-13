/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.util.Random;

/**
 * Multivalent is the abstract base class for Bivalent and Quadrivalent
 * Note that in all constructors of descendant classes popdata and haplostruct
 * must be set!
 * @author Roeland Voorrips
 */
abstract public class Multivalent {

    public static final double RECOMBDIST = 0.5;

    //3 protected fields for faster access, must be set by constructor
    //of descendant classes:
    protected PopulationData popdata;
    protected Tools tools;
    protected Random rand;

    protected Chromosome chrom;

    /**
     * haplostruct are the two (Bivalent) or 4 (Quadrivalent) chromosomes
     * that make up the multivalent. Their order is NOT randomized in the
     * constructors
     */
    protected HaploStruct[] haplostruct;
        // the 2 or 4 chromosomes involved in a Bivalent or Quadrivalent
        // (before recombination occurs)

    public HaploStruct[] getHaplostruct() {
        return haplostruct;
    }

    public Chromosome getChrom() {
        return chrom;
    }

    /**
     * doMeiosis performs a meiosis. This is a simplified
     * model in that chromatid interference is not modelled.
     * The recombinations take place independently of each other, each between a
     * random pair of chromatids belonging to different chromosomes.
     * In the first meiotic division half of the centromeres are assigned
     * at random to the first and second pair of gametes, meaning that
     * gamete[0] is a random gamete. Gamete pairs 0&1, and 2&3, have the same
     * centromere allele(s).
     *
     * doMeiosis is a three-step process that internally calls:
     * 1 - doCrossingOver (produces 2 recombined "slots" for each haplostruct)
     * 2 - slotsToPatterns (produces one "pattern" for each slot)
     * 3 - patternsToGametes (produces four gametes from the patterns)
     * A slot represents a fixed position in the multivalent containing one
     * of the chromatids. It uses a HaploStruct to record at which positions
     * it connects to which other slots (the founder component of the
     * slot HaploStruct contains the target slots).
     * A pattern represents the mosaic of slot segments that results after
     * the crossing-over process is finished. A pattern is also implemented
     * as a HaploStruct, where the founder component records the slot
     * numbers.
     * Each pattern is finally converted to a chromatid in a gamete by
     * applying the pattern of slot segments to the original composition
     * of the chromosomes that filled those slots.
     * @return a HaploStruct[] with 2 elements for each chromosome in the
     * Multivalent: the recombined chromatids (4 for a Bivalent, 8 for a
     * Quadrivalent ordered per gamete: 0+1 in gamete 0, 2+3 in gamete 1 etc.)
     * The gametes are ordered as a tetrad, with the first two deriving
     * from the "top" side of the first meiotic division and the last two
     * from the "bottom" side. Which chromosomes (centromeres) make up
     * the top and bottom sides is determined randomly.
     * @throws Exception
     */
    public HaploStruct[] doMeiosis() throws Exception {
        return patternsToGametes(slotsToPatterns(doCrossingOver()),
                getCentromereSortOrder());
    } //doMeiosis

    public double distToFirstRecomb (double meandist) {
        double pos = tools.ranExp(meandist);
        if (popdata.chiasmaInterference) {
            double burninDistance = 10.0*meandist; //to obtain stationary distribution, burn-in 10 times
            while (pos<burninDistance) {
                pos += tools.ranDistInterference(meandist);
            }
            return pos-burninDistance;
        }
        else {
            // no interference
            return pos;
        }
    } //distToFirstRecomb

    public double distToNextRecomb (double meandist) {
        if (popdata.chiasmaInterference) {
            return tools.ranDistInterference(meandist);
        } else {
            return tools.ranExp(meandist);
        }
    } //distToNextRecomb

    /**
     * doCrossingOver :
     * assume the (2*ploidy) chromatids to be fixed in parallel "slots",
     * so that slot 0+1 contain the chromosome 0 over their full length,
     * slot 2+3 have chromosome 1 (and so on for polyploids). 
     * A recombination occurs between any two non-sister slots;
     * we record its position and the involved slots, but the content of each
     * slot over its full length remains unchanged. In this implementation
     * intra-chromosomal recombinations (between sister slots) are always without
     * effect and therefore can be ignored (and we don't generate them).
     * We (mis-)use the HaploStruct class to record the recombinations for
     * each slot: each "segment" now has the recombination position and the
     * slot it connects to.
     * (NOTE that this is realistic for a bivalent but not for a
     * quadrivalent; in Quadrivalent this method is overridden).
     * In slotsToPatterns this information is used to produce "patterns":
     * HaploStructs that describe how each resulting chromatid is composed
     * from fragments in the different slots.
     * In patternsToGametes these patterns are separated according to their
     * centromeric composition to form the four gametes
     * (In the first meiotic division the first half and the second half
     * of the patterns are separated from each other, and in the second division
     * one pattern of each subsequent pair is moved to one of the two gametes.)
     * Finally (in fillGamete) the patterns in the gametes are used to produce
     * the actual recombined chromatids, by applying the patterns to the
     * original chromosomes which were themselves HaploStructs composed of
     * founder alleles.
     *
     * @return the 2*ploidy slots after recombination, still in the same order
     * (i.e. slots 0+1 for the first chromosome etc)
     * @throws Exception 
     */
    protected HaploStruct[] doCrossingOver() throws Exception {
        //update counts (only used in test runs):
        if (popdata.ploidy > 2) {
            if (haplostruct.length==2) ((TetraploidChromosome)chrom).bivalentCount++;
            else ((TetraploidChromosome)chrom).paralQuadrivalentCount++;
        } 
        
        if (popdata.bivalentsBidirectional) {
            return doBidirectionalCrossingOver();
        }
        //else normal CrossingOver:
        int slotcount = haplostruct.length * 2; // 2 chromatids for each chromosome
        double recombDist = (haplostruct.length==2) ?
            RECOMBDIST : RECOMBDIST * 2 / haplostruct.length;
            // i.e. more (parallel) chromosomes require shorter 
            // recombination distance
        if (!popdata.allowNoRecomb) {
            recombDist = recombDist/0.5 * Tools.calcRecombDist(chrom.getLength());
        }
        HaploStruct[] slots = new HaploStruct[slotcount];
        do { //repeat once or (if allowNoRecomb false) until at least one recombination
            for (int s=0; s<slotcount; s++) {
                slots[s] = new HaploStruct(chrom, -1); //the first segment will be ignored
            }
            double pos =  //position of first recombination
                chrom.getHeadPos() + distToFirstRecomb(recombDist);
            int[] ix = new int[2]; //the two slots for the currect recombination
            while (pos < chrom.getTailPos()) {
                //loop: apply crossing-over and get next recombination
                ix[0] = rand.nextInt(slotcount);
                do {
                    ix[1] = rand.nextInt(slotcount);
                } while ( ix[1]/2 == ix[0]/2 ); //no intra-chromosomal recombinations
                /* ix[0] and ix[1] are now the chromatids involved in crossing-over,
                 * to each we add a "segment" with recombination position and the slot
                 * it connects to:
                 */
                slots[ix[0]].addSegment(pos,ix[1]);
                slots[ix[1]].addSegment(pos,ix[0]);
                pos += distToNextRecomb(recombDist);
            } //while pos
        } while (!(popdata.allowNoRecomb || allChromRecomb(slots)));
        return slots;
    } //doCrossingOver

    /**
     * rejectDist:
     * in a cross-type Quadrivalent recombinations are generated from both ends
     * of each chromosome. A new recombination generated from one end may not "fit"
     * as seen from the other end: it may be located beyond that end of the
     * chromosome, or it may be located beyond the last recombination generated
     * from that end. If interference occurs it may also be located too close
     * to the last recombination from that end.
     * We use this function rejectDist to tell whether the new recombination position
     * is acceptable, depending on the popdata.chiasmaInterference setting.
     * In order to validate this function we tested it in bivalents, where
     * we compared the results of generating the recombinations from one end
     * with those obtained generating recombinations from both ends (by setting
     * popdata.bivalentsBidirectional to true and so invoking
     * Multivalents.doBidirectionalCrossingOver).
     * Experimentation showed that the following procedure gives results
     * indistinguihable from the unidirectionally generated recombinations
     * (i.e. both the means and the standard deviations of the obtained
     * recombination frequencies are very close to the unidirectional ones):
     * - if popdata.chiamaInterference is false, the new recombination is
     *   rejected if it lies beyond the last one generated from the other side.
     * - if chiasmaInterference is true and the distance from the new recombination
     *   to the last recombination from the opposite side is less than 0.5 Morgan,
     *   the probability of rejection increases from 0 at a distance of 0.5 M 
     *   to 1 at 0 M, following a slightly steeper than quadratic function
     *   (power = 2.25).
     * @param dist the distance remaining between the new recombination and
     * the last recombination generated from the other end; negative if the 
     * new recombination lies beyond (at the wrong side) of the existing 
     * recombination
     * @return true if rejected
     */
    protected boolean rejectDist(double dist) {
        return (dist<=0 ||
               (popdata.chiasmaInterference && (dist < 0.5) && 
                   (rand.nextDouble() < Math.pow(1.0-2*dist, 2.25))));
    } //rejectDist

    /**
     * doBidirectionalCrossingOver does the same as doCrossingOver,
     * but the recombinations are generated alternating from the head and
     * the tail of the chromosome.
     * This is the way they are generated in Quadrivalents (starting from
     * the ends of each arm) because that is the only logical way in that case,
     * but we run into the decision whether a new recombination can or cannot be
     * added between two existing ones, especially in the case of chiasma
     * interference.
     * This procedure is added to test if in a bivalent both methods lead
     * to the same result. (They do indeed!)
     * @return as in doCrossingOver
     * @throws Exception
     */
    private HaploStruct[] doBidirectionalCrossingOver() throws Exception {
        int slotcount = haplostruct.length * 2;
        double recombDist = (haplostruct.length==2) ?
            RECOMBDIST : RECOMBDIST * 2 / haplostruct.length;
            //i.e. more (parallel) chromosomes require shorter recombination distance
        if (!popdata.allowNoRecomb) {
            recombDist = recombDist/0.5 * Tools.calcRecombDist(chrom.getLength());
        }
        HaploStruct[] slots = new HaploStruct[slotcount];
        do { //repeat once or (if allowNoRecomb false) until at least one recombination)
            for (int s=0; s<slotcount; s++) {
                slots[s] = new HaploStruct(chrom, -1); //the first segment will be ignored
            }
            double[] pos = new double[] {Double.NaN, Double.NaN};
               //new position of next recombination from chromosome head and tail
            int[] ix = new int[2]; //the two slots for the currect recombination

            /* First recombination, from head or from tail, outside loop
             * because it is tested against the other end, not against another
             * recombination, and therefore interference is not an issue
             */
            int side = rand.nextInt(2); //0=head of chrom, 1 = tail of chrom
            boolean finished = false; //true when no more recombinations can be added
            if (side==0) {
                //start from chromosome head
                pos[0] = chrom.getHeadPos() + distToFirstRecomb(recombDist);
                if (pos[0] >= chrom.getTailPos()) finished=true;
            }
            else { //side=1
                //start from chromosome tail
                pos[1] = chrom.getTailPos() - distToFirstRecomb(recombDist);
                if (pos[1] <= chrom.getHeadPos()) finished=true;
            }
            if (!finished) {
                //add the first recombination
                ix[0] = rand.nextInt(slotcount);
                do {
                    ix[1] = rand.nextInt(slotcount);
                } while ( ix[1]/2 == ix[0]/2 ); //no intra-chromosomal recombinations
                slots[ix[0]].addSegment(pos[side],ix[1]);
                slots[ix[1]].addSegment(pos[side],ix[0]);
            }

            // all further recombinations are added in a loop:
            while (!finished) {
                side = 1-side;
                if (Double.isNaN(pos[side])) { //occurs only at the second recombination
                    if (side==0) {
                        pos[0] = chrom.getHeadPos() + distToFirstRecomb(recombDist);
                    }
                    else {
                        pos[1] = chrom.getTailPos() - distToFirstRecomb(recombDist);
                    }
                }
                else {
                    pos[side] = pos[side] - (2 * side - 1) * distToNextRecomb(recombDist);
                }
                if (rejectDist((2*side-1)*(pos[side]-pos[1-side]))) finished=true;
                else {
                    ix[0] = rand.nextInt(slotcount);
                    do {
                        ix[1] = rand.nextInt(slotcount);
                    } while ( ix[1]/2 == ix[0]/2 ); //no intra-chromosomal recombinations
                    slots[ix[0]].insertSegment(pos[side],ix[1]);
                    slots[ix[1]].insertSegment(pos[side],ix[0]);
                }
            }
        } while (!(popdata.allowNoRecomb || allChromRecomb(slots)));
        return slots;
    } //doBidirectionalCrossingOver



    /**
     * slotsToPatterns :
     * In slots, the founder array contains references to which other
     * slots each recombination connects (see doCrossingOver).
     * As each chromosome consists of two chromatids that are positioned
     * in two slots, in a bivalent these "founders" can range from 0 to 3
     * (4 slots) and in a quadrivalent from 0 to 7 (8 slots).
     * In slotsToPatterns, the branching patterns stored in the slots is
     * converted into a set of mosaics ("patterns") that are obtained
     * by following each recombined chromatid from the head of a slot
     * through the branching pattern. Therefore the patterns are composed
     * of a sequence of the slots from which they are made up. But the slot
     * numbers are converted into the chromosome number they belong to
     * (in a bivalent: slots 0,1,2,3 -> chromosome 0,0,1,1, similar for a
     * quadrivalent). Therefore the founder array of a pattern can contain
     * only numbers 0 and 1 (bivalent) or 0..3 (quadrivalent).
     * Note that the order of the patterns is determined by the chromosome
     * head: patterns 0 and 1 start with founder 0, patterns 2 and 3 with
     * founder 1 etc.
     * @param slots
     * @return
     * @throws Exception
     */
    private HaploStruct[] slotsToPatterns(HaploStruct[] slots) throws Exception {
        int slotcount = slots.length;
        //for (int s=0; s<slotcount; s++) System.out.println("slots["+s+"]:\t"+slots[s]);
        HaploStruct[] patterns = new HaploStruct[slotcount];
        double pos;
        for (int p=0; p<slotcount; p++) {
            patterns[p] = new HaploStruct(chrom, p/2); //first segment has "founder" p/2: 0,1,2,3 -> 0,0,1,1
            int currslot = p; //the current slot;
            int nextseg = 1; //the next segment on the current slot
            while (nextseg < slots[currslot].segmentCount()) {
                pos = slots[currslot].getRecombPos().get(nextseg);
                int slt = slots[currslot].getFounder().get(nextseg); //next slot
                /* add the chromosome corresponding to slot slt:
                 * 0+1->0, 2+3->1, etc for polyploids
                 */
                patterns[p].addSegment(pos, slt/2);
                currslot = slt;
                nextseg = slots[currslot].findSegment(pos);
                if (slots[currslot].getRecombPos().get(nextseg) != pos) {
                    nextseg = Integer.MAX_VALUE;
                }
                nextseg++;
            }
        } //for pb
        return patterns;
    } //slotsToPatterns
    
    /** 
     * getCentromereSortOrder produces a list of centromere alleles
     * (0 .. haplostruct.length-1, each occurring twice) sorted such that
     * they reflect the gamete order (see patternsToGametes).
     * In short:
     * the centromeres are divided in two groups at the first meiotic division;
     * they are then each split in half and these two halves are separated
     * in the second meiotic division.
     * The first two gametes therefore have the same centromere alleles,
     * and the second two gametes likewise.
     * Note: As the chromosomes are NOT randomized in the Bi- and Quadrivalent 
     * constructors, the random assignment of centromeres to the two poles
     * of the first meiotic division must be done here. 
     * Question 15-4-13: is that still true? Randomization not in constructor
     * but in Individual.calcChromConfig?
     * Answer: the 2 or 4 chromosomes are randomized already
     * in Individual.calcChromConfig, but for cross-type Quadrivalents
     * separation in first meiotic division should be independent of the 
     * arms structure;
     * therefore randomization of centromeres still seems important
     * (and even if it is not, it is certainly not wrong)
     * The order of the chromosomes in each gamete doesn't really matter
     * but to make it as random as possible the order of the chromosomes
     * in each gamete is randomized.
     * This procedure guarantees that:
     * - The four gametes are nicely ordered in a tetrad, with the first two
     *   and the second two gametes derived from opposite poles in the 
     *   first meiotic division
     * - The order of the homologous chromosomes if completely independent
     *   between any pair of gametes
     * - the first gamete is a random gamete: taking always the first gamete
     *   from successive meioses of the same parent will yield a random
     *   sample of gametes
     * 
     * @return an int array specifying the order of the centromere alleles
     * with 2* the length of the haplostruct array; each quarter of the
     * array representing the chromosomes in one gamete.
     */
    int[] getCentromereSortOrder() {
        int[] centro;
        if (haplostruct.length==2) {
            //bivalent
            if (rand.nextBoolean()) 
                 centro = new int[] {0,0,1,1};
            else centro = new int[] {1,1,0,0};
        }
        else {
            //parallel multivalent (Quadrivalent overrides this method)
            int[] chr = tools.shuffleArray(Tools.seq(haplostruct.length));
            /* for a 4-valent: {a,b,c,d}, 6-valent: {a,b,c,d,e,f) etc
             * with a:d = randomized 0:3 or a:f = randomized 0:5
             * The first meiotic division separates ab/cd, abc/def etc:
             */
            int[][] gam = new int[4][haplostruct.length/2];
            /* gam[0..4] are the 4 gametes and will contain the ORDERED
             * sets of centromere alleles
             */
            for (int p=0; p<haplostruct.length/2; p++) {
                gam[0][p] = chr[p];
                gam[1][p] = chr[p];
                gam[2][p] = chr[p+haplostruct.length/2];
                gam[3][p] = chr[p+haplostruct.length/2];
            }
            /* we have now for 4-valent: { {a,b}, {a,b}, {c,d}, {c,d} },
             * for 6-valent: { {a,b,c}, {a,b,c}, {d,e,f}, {d,e,f} }
             * and similar for higher multivalents
             */
            // randomize the order in each gamete:
            for (int g=0; g<4; g++) {
                gam[g] = tools.shuffleArray(gam[g]);
            }
            //next, put gam[0..3] consecutively into centro:
            centro = new int[2*haplostruct.length];
            for (int g=0; g<4; g++) {
                System.arraycopy(gam[g], 0, centro, g*(centro.length/4), centro.length/4);
            }
        } //quadrivalent or higher
        return centro;
    } //getCentromereSortOrder

    /**
     * patternsToGametes :
     * The patterns array stores a description of how the recombinants
     * are composed of sections of the original chromosomes of the
     * Multivalent. Two more actions need to be performed to convert
     * the patterns array to a gametes array:
     * (1)the patterns must be ordered such that the gametes array
     *    reflects the order of the two meiotic divisions:
     *    - the first meiotic division separates entire centromeres,
     *      one half of the centromeres (1 in a bivalent, 2 in a quadrivalent)
     *      going to each pole;
     *    - the second meiotic division splits each centromere in two halves
     *      and assigns each half at random to one of the two gametes at that 
     *      pole.
     *    - the total of 4 (bivalent) or 8 (quadrivalent) haplotypes in the
     *      gametes array are ordered in the order of the gametes, i.e.
     *      * for a bivalent (4 haplotypes) nr 0, 1, 2 and 3 each go a different 
     *        gamete, where 0 and 1 represent one pole of the first division
     *        and 2 and 3 the other pole;
     *      * for a quafrivalent, nrs 0 and 1 are in the first gamete, 
     *        2 and 3 in the second etc; where again the first two gametes
     *        are from one pole and the last two gametes from the other pole.
     *    - the sorting of the patterns haplotypes array therefore takes place
     *      based of the allele in founder array at the centromere position.
     *      Question 15-4-13: Is that correct? except in a founder several homologs
     *      may have the same founder allele at the centromere ?
     *      Answer: although it is called the founder array, in a pattern
     *      haplotype the founder array indicates the slots, not the real
     *      founder allele. So yes, it is correct.
     * (2) while the patterns describe how the recombinants are composed of
     *     sections of the original 2 (bivalent) or 4 (quadrivalent) chromosomes
     *     in the multivalent, the gametes must be specified in terms of the
     *     original founder alleles of the population; i.e. the composition of
     *     the multivalent chromosomes must be taken into account.
     * @param patterns a HaploStruct array of length 2*ploidy describing the
     *     recombinants as mosaics of the original chromosomes of the multivalent.
     * @param centro is an array of 4 (bivalent) or 8 (quadrivalent) numbers
     *     that specifies the order of the centromere alleles in the
     *     gametes array
     * @return a similar HaploStruct array of length 2*ploidy containing the
     * actual chromatids in the gametes; sorted such that the (1 or 2)
     * chromatids of each gamete are consecutive, and that the first half
     * of the chromatids (the first two gametes) and the second half (the
     * second two gametes) reflect the separation in the first meiotic
     * division.
     * @throws Exception
     */
    private HaploStruct[] patternsToGametes(HaploStruct[] patterns, int[] centro)
            throws Exception {
        //sort the patterns such that the four gametes
        //are based on the centromere alleles, in the order of centro:
        HaploStruct tmphs;
        for (int s=0; s<patterns.length-1; s++) {
            if (patterns[s].getFounderAtCentromere() != centro[s]) {
                int t = s+1;
                while (patterns[t].getFounderAtCentromere() != centro[s]) {
                    t++;
                }
                tmphs = patterns[s];
                patterns[s] = patterns[t];
                patterns[t] = tmphs;
            }
        }
        /*
         * Note that in a quadrivalent the centro array contains two 0's, 
         * two 1's, two 2's and two 3's, and that the two chromatids with
         * the same centromere allele are not yet randomized. We do that now:
         */
        if (patterns.length>4) {
            for (int all=0; all<4; all++) {
                int first, second;
                if (rand.nextBoolean()) {
                    first=0; while (centro[first]!=all) first++;
                    second=first+1; while (centro[second]!=all) second++;
                    tmphs = patterns[first];
                    patterns[first] = patterns[second];
                    patterns[second] = tmphs;
                }
            }
        }
        

        //finally construct the gametes based on the slots structure:
        HaploStruct[] gametes = new HaploStruct[patterns.length];
        for (int s=0; s<patterns.length; s++) {
            gametes[s] = fillGamete(patterns[s]);
        }
        return gametes;
    } //patternsToGametes

    /**
     * fillGamete takes a pattern: a HaploStruct composed of 0's and 1's
     * (and higher numbers if this Multivalent is not a Bivalent)
     * that refer to this Multivalents's haplostruct
     * and uses that to fill a gamete's HaploStruct in terms of
     * founder alleles.
     * @param pattern
     * @return
     * @throws Exception
     */
    private HaploStruct fillGamete(HaploStruct pattern) throws Exception {
        HaploStruct gamete;
        if (pattern.segmentCount()==1) {
            gamete = (HaploStruct) haplostruct[pattern.getFounder().get(0)].clone();
        }
        else {
            gamete = new HaploStruct(getChrom(),pattern.getFounder().get(0));
            gamete.getFounder().clear();
            gamete.getRecombPos().clear();
            int patseg = 0; //the pattern segment currently being translated
            while (patseg<pattern.segmentCount()) {
                /* get the list of partial and complete segments
                 * of the haplostruct corresponding to this pattern segment;
                 * the start of the first (partial) segment is replaced
                 * by the start of the pattern segment
                 * and all the haplostruct segments are added to gamete
                 */
                double patstart = pattern.getRecombPos().get(patseg);
                double patend = patseg < pattern.segmentCount()-1 ?
                    pattern.getRecombPos().get(patseg+1) :
                    getChrom().getTailPos();
                int pf = pattern.getFounder().get(patseg); //always in 0..3 for Quadrivalent, 0..1 for Bivalent (?)
                HaploStruct hs = haplostruct[pf];
                //find the hsseg in which patstart falls:
                int hsseg = hs.findSegment(patstart);
                gamete.getRecombPos().add(patstart);
                gamete.getFounder().add(hs.getFounder().get(hsseg));
                hsseg++;
                while (hsseg<hs.segmentCount() &&
                        hs.getRecombPos().get(hsseg)<patend) {
                    gamete.getRecombPos().add(hs.getRecombPos().get(hsseg));
                    gamete.getFounder().add(hs.getFounder().get(hsseg));
                    hsseg++;
                }
                patseg++;
            }
        }
        return gamete;
    } //fillGamete
    
    /**
     * allChromRecomb:
        checks if at least one chromatid of each chromosome is involved in
        at least one recombination
     * @param slots each chromatid (a HaploStruct) occupies one slot;
     * a Bivalent has 4 slots, a Quadrivalent 8 slots
     * @return true if the check is positive
     */
    boolean allChromRecomb(HaploStruct[] slots) {
        if (slots.length==4) { 
            //bivalent: only check one chrom, if this has recomb then the other has too
            return slots[0].getFounder().size()>1 ||
                   slots[1].getFounder().size()>1;
        }
        else {
            boolean result;
            int p = slots.length/2 - 1;
            do {
                result = slots[2*p].getFounder().size()>1 ||
                         slots[2*p+1].getFounder().size()>1;
                p--;
            } while (result && p>=0);
            return result;
        }
    } //allChromRecomb
    
}

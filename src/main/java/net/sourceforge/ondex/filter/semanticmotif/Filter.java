package net.sourceforge.ondex.filter.semanticmotif;

import net.sourceforge.ondex.args.ArgumentDefinition;
import net.sourceforge.ondex.filter.ArgumentNames;
import net.sourceforge.ondex.args.BooleanArgumentDefinition;
import net.sourceforge.ondex.args.StringArgumentDefinition;
import net.sourceforge.ondex.core.*;
import net.sourceforge.ondex.core.util.DefaultBitSet;
import net.sourceforge.ondex.core.util.ONDEXBitSet;
import net.sourceforge.ondex.core.util.ONDEXViewFunctions;
import net.sourceforge.ondex.core.util.SparseBitSet;
import net.sourceforge.ondex.event.type.GeneralOutputEvent;
import net.sourceforge.ondex.event.type.WrongParameterEvent;
import net.sourceforge.ondex.filter.ONDEXFilter;
import net.sourceforge.ondex.tools.ondex.ONDEXGraphCloner;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Map;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.ArrayList;

/**
 * Finds Semantic motifs 
 *
 * @author taubertj
 * @version 31.01.2008
 */
public class Filter extends ONDEXFilter implements ArgumentNames {

    // contains list of visible concepts
    private Set<ONDEXConcept> concepts = null;

    // contains list of visible relations
    private Set<ONDEXRelation> relations = null;

    // contains list of invisible concepts
    private Set<ONDEXConcept> invconcepts = null;

    // contains list of invisible relations
    private Set<ONDEXRelation> invrelations = null;

    /**
     * Constructor
     */
    public Filter() {
    }

    @Override
    public void copyResultsToNewGraph(ONDEXGraph exportGraph) {
        ONDEXGraphCloner graphCloner = new ONDEXGraphCloner(graph, exportGraph);
        for (ONDEXConcept c : concepts) {
            graphCloner.cloneConcept(c);
        }
        for (ONDEXRelation r : relations) {
            graphCloner.cloneRelation(r);
        }
    }

    @Override
    public Set<ONDEXConcept> getVisibleConcepts() {
        return concepts;
    }

    @Override
    public Set<ONDEXRelation> getVisibleRelations() {
        return relations;
    }

    public Set<ONDEXConcept> getInVisibleConcepts() {
        return invconcepts;
    }

    public Set<ONDEXRelation> getInVisibleRelations() {
        return invrelations;
    }

    /**
     * Returns the name of this filter.
     *
     * @return name
     */
    public String getName() {
        return "SemanticMotif Filter";
    }

    /**
     * Returns the version of this filter.
     *
     * @return version
     */
    public String getVersion() {
        return "20.01.2010";
    }

    @Override
    public String getId() {
        return "semanticmotif";
    }


    /**
     * No Arguments (yet)
     *
     * @return single argument definition
     */
    public ArgumentDefinition<?>[] getArgumentDefinitions() {
        return new ArgumentDefinition<?>[0];
    }

    public String[] getBindingPIDs(Set<ONDEXRelation> drugRelations, ONDEXConcept drugConcept) {
        ArrayList pids = new ArrayList();
        for (ONDEXRelation dr : drugRelations) {
            if (dr.getOfType().getId().equals("bi_to")) {
                ONDEXConcept from = dr.getFromConcept();
                ONDEXConcept to = dr.getToConcept();
                Set<ConceptAccession> it;
                if (from.equals(drugConcept)) {
                    //pids.add(to.getConceptName().getName());
                    //pids.add(to.getPID());
                    it = to.getConceptAccessions();
                }
                else {
                    //pids.add(from.getConceptName().getName());
                    //pids.add(from.getPID());
                    it = from.getConceptAccessions();
                }
                for (ConceptAccession acc : it) {
                    if (acc.getElementOf().equals(graph.getMetaData().getCV("uniprot"))) {
                        pids.add(acc.getAccession());
                    }
                }
            }
        }
        String[] pid_array = new String[pids.size()];
        for (int i = 0; i < pids.size(); i++) {
            pid_array[i] = (String) pids.get(i);
        }
        return pid_array;
    }

    //public boolean checkOverlap(String[] b1, String[] b2) {
    public int checkOverlap(String[] b1, String[] b2) {
        int i = 0;
        boolean match = true;
        for (String id1 : b1) {
            boolean temp_match = false;
            for (String id2 : b2) {
                System.out.println(id1+"\t"+id2);
                if (id1.equals(id2)) {
                    temp_match = true;
                }
            }
            if (temp_match == false) {
                i++;
                match = false;
            }
        }
        System.out.println(i+"");
        /*if (match == false && b1.length != 0 && b2.length != 0) {
            return true;
        }
        else {
            return false;
        }*/
        return i;
    }

    /**
     * Filters the graph and constructs the lists for visible concepts and
     * relations.
     */
    public void start() {
        int motifs = 0;
        concepts = ONDEXViewFunctions.copy(graph.getConcepts());
        relations = ONDEXViewFunctions.copy(graph.getRelations());
        Set<ONDEXRelation> relationview = graph.getRelationsOfRelationType(graph.getMetaData().getRelationType("sim"));
        for (ONDEXRelation r : relationview) {
            //get the 2 drugs connected to this relation
            ONDEXConcept drug1 = r.getFromConcept();
            ONDEXConcept drug2 = r.getToConcept();
            Set<ONDEXRelation> drug1Relations = graph.getRelationsOfConcept(drug1);
            Set<ONDEXRelation> drug2Relations = graph.getRelationsOfConcept(drug2);
            //find ids of target concepts, store in array
            String[] bind1 = getBindingPIDs(drug1Relations, drug1);
            String[] bind2 = getBindingPIDs(drug2Relations, drug2);
            //check for overlap
            motifs = motifs + checkOverlap(bind1, bind2);
            /*boolean motif = checkOverlap(bind1, bind2);
            if (motif == true) {
                motifs++;
            }*/
        }
        System.out.println(motifs);
        /*
        ONDEXView<ONDEXConcept> conceptview = graph.getConceptsOfConceptClass(graph.getMetaData().getConceptClass("Compound"));
        while(conceptview.hasNext()) {
            ONDEXConcept c1 = conceptview.next();
            //list of similar drugs ONDEXConcept c2 = null;
            Map<String, ONDEXConcept> sims = new Hashtable<String, ONDEXConcept>();
            //list of targets
            Map<String, ONDEXConcept> targets = new Hashtable<String, ONDEXConcept>();
            ONDEXView<ONDEXRelation> relationsofc1 = graph.getRelationsOfConcept(c1);
            while (relationsofc1.hasNext()) {
                ONDEXRelation r = relationsofc1.next();
                ONDEXConcept from = r.getFromConcept();
                ONDEXConcept to = r.getToConcept();
                RelationType rt = r.getOfType();
                if (rt.getId().equals("sim")) {
                    if (from.equals(c1)) {
                        sims.put(to.getPID(), to);
                    }
                    else {
                        sims.put(from.getPID(), from);
                    }
                }
                else if (rt.getId().equals("bi_to")) {
                    if (from.equals(c1)) {
                        targets.put(to.getPID(), to);
                    }
                    else {
                        targets.put(from.getPID(), from);
                    }
                }
            }
            Iterator e = sims.keySet().iterator();
            while (e.hasNext()) {
                boolean match = false;
                ONDEXConcept c2 = sims.get(e.next());
                Map<String, ONDEXConcept> sim_targets = new Hashtable<String, ONDEXConcept>();
                ONDEXView<ONDEXRelation> relationsofc2 = graph.getRelationsOfConcept(sims.get(c2));
                while (relationsofc2.hasNext()) {
                    ONDEXRelation r = relationsofc2.next();
                    ONDEXConcept to = r.getToConcept();
                    ONDEXConcept from = r.getFromConcept();
                    RelationType rt = r.getOfType();
                    if (rt.getId().equals("bi_to")) {
                        if (from.equals(c2)) {
                            sim_targets.put(to.getPID(), to);
                        }
                        else {
                            sim_targets.put(from.getPID(), from);
                        }
                    }
                }
                //does this similar drug share any targets with the query drug?
                Iterator e1 = targets.keySet().iterator();
                while (e1.hasNext()) {
                    Iterator e2 = sim_targets.keySet().iterator();
                    String key1 = e1.next().toString();
                    while (e2.hasNext()) {
                        String key2 = e2.next().toString();
                        if (key1.equals(key2)) {
                            match = true;
                        }
                    }
                }
                if (match == false) {
                    System.out.println("Semantic Motif! " + c1.getConceptName() + "    " + c2.getConceptName());
                }
            }
        }
        */
    }

    /**
     * An indexed graph is not required.
     *
     * @return false
     */
    public boolean requiresIndexedGraph() {
        return false;
	}

	public String[] requiresValidators() {
		return new String[0];
	}
}

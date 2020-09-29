# strings for sample queries

function path_pc(pathVals::String)
    """
    PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
    PREFIX dcterms: <http://purl.org/dc/terms/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>

    SELECT DISTINCT ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller (GROUP_CONCAT(?leref; SEPARATOR = ",") AS ?lerefs) (GROUP_CONCAT(?reref; SEPARATOR = ",") AS ?rerefs)

    WHERE {
        OPTIONAL { ?left a bp:Complex ;
                         bp:xref [ a bp:UnificationXref ;
                                   bp:db "reactome"^^xsd:string ;
                                   bp:id ?leref] }

        OPTIONAL { ?right a bp:Complex ;
                         bp:xref [ a bp:UnificationXref ;
                                   bp:db "reactome"^^xsd:string ;
                                   bp:id ?reref] }

        OPTIONAL { ?left a bp:Protein ;
                         bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
                                                                               bp:db ?dbname ;
                                                                               bp:id ?leref ]
        			FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?right a bp:Protein ;
                          bp:entityReference/bp:memberEntityReference/bp:xref [ a bp:UnificationXref ;
                                                                                bp:db ?dbname ;
                                                                                bp:id ?reref ]
        			FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?left a bp:Protein ;
                         bp:entityReference [ a bp:ProteinReference ;
                                              bp:xref [ bp:db ?dbname ;
                                                        bp:id ?leref ] ] .
                    FILTER ( NOT EXISTS {?left bp:memberPhysicalEntity ?lmem } )
                    FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?right a bp:Protein ;
                          bp:entityReference [ a bp:ProteinReference ;
                                               bp:xref [ bp:db ?dbname ;
                                                         bp:ed ?reref ] ] .
                    FILTER ( NOT EXISTS {?right bp:memberPhysicalEntity ?rmem } )
                    FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?left a bp:Dna ;
                          bp:entityReference [ a bp:DnaReference ;
                                               bp:xref [ bp:db ?dbname ;
                                                         bp:id ?leref] ]
    				FILTER regex(?dbname, "^.*uniprot") }

    	OPTIONAL { ?right a bp:Dna ;
                          bp:entityReference [ a bp:DnaReference ;
                                               bp:xref [ bp:db ?dbname ;
                                                         bp:id ?leref] ]
    				FILTER regex(?dbname, "^.*uniprot") }

        {
            SELECT ?pw ?pathname ?pc ?ctype
            		?left ?right ?ltype ?rtype ?llocref ?rlocref
            		?cltype ?controller
            WHERE {
                ?pw bp:name ?pathname ;
                    bp:pathwayComponent* ?pc .
            VALUES ?pw $pathVals
                ?pc a ?ctype .
                ?pc a bp:Interaction .
                ?left ^bp:left ?pc ;
                       a ?ltype ;
                       bp:cellularLocation [ bp:term ?lloc ;
                                             bp:xref ?llocref ] ;
                       bp:displayName ?lname .
                ?right ^bp:right ?pc ;
                       a ?rtype ;
                       bp:cellularLocation [ bp:term ?rloc ;
                                             bp:xref ?rlocref ] ;
                       bp:displayName ?rname .
                OPTIONAL { ?pc ^bp:controlled [ bp:controlType ?cltype ;
                                                bp:controller ?controller ]}
            }
        }
        VALUES ?rtype {bp:Dna bp:Protein bp:Complex}
        VALUES ?ltype {bp:Dna bp:Protein bp:Complex}
        VALUES ?ctype {bp:BiochemicalReaction}
    }

    GROUP BY ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller
    """
end

function path_pc_old(pathVals::String)
    """
    PREFIX wp:    <http://vocabularies.wikipathways.org/wp#>
    PREFIX dcterms: <http://purl.org/dc/terms/>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>

    SELECT DISTINCT ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller (GROUP_CONCAT(?leref; SEPARATOR = ",") AS ?lerefs) (GROUP_CONCAT(?reref; SEPARATOR = ",") AS ?rerefs)

    WHERE {
        OPTIONAL { ?left a bp:Protein ;
                         bp:entityReference [ a bp:ProteinReference ;
                                              bp:xref [ bp:db ?dbname ;
                                                        bp:id ?leref ] ] .
                    FILTER ( NOT EXISTS {?left bp:memberPhysicalEntity ?lmem } )
                    FILTER regex(?dbname, "^.*uniprot") }

        OPTIONAL { ?right a bp:Protein ;
                          bp:entityReference [ a bp:ProteinReference ;
                                               bp:xref [ bp:db ?dbname ;
                                                         bp:ed ?reref ] ] .
                    FILTER ( NOT EXISTS {?right bp:memberPhysicalEntity ?rmem } )
                    FILTER regex(?dbname, "^.*uniprot") }

        {
            SELECT ?pw ?pathname ?pc ?ctype
            		?left ?right ?ltype ?rtype ?llocref ?rlocref
            		?cltype ?controller
            WHERE {
                ?pw bp:name ?pathname ;
                    bp:pathwayComponent* ?pc .
            VALUES ?pw $pathVals
                ?pc a ?ctype .
                ?pc a bp:Interaction .
                ?left ^bp:left ?pc ;
                       a ?ltype ;
                       bp:cellularLocation [ bp:term ?lloc ;
                                             bp:xref ?llocref ] ;
                       bp:displayName ?lname .
                ?right ^bp:right ?pc ;
                       a ?rtype ;
                       bp:cellularLocation [ bp:term ?rloc ;
                                             bp:xref ?rlocref ] ;
                       bp:displayName ?rname .
                OPTIONAL { ?pc ^bp:controlled [ bp:controlType ?cltype ;
                                                bp:controller ?controller ]}
            }
        }
        VALUES ?rtype {bp:Dna bp:Protein bp:Complex}
        VALUES ?ltype {bp:Dna bp:Protein bp:Complex}
        VALUES ?ctype {bp:BiochemicalReaction}
    }

    GROUP BY ?left ?right ?pw ?pathname ?pc ?ctype ?ltype ?rtype ?llocref ?rlocref ?cltype ?controller
    """
end

function go_obo(terms::String)
    """
    PREFIX obo-term: <http://purl.obolibrary.org/obo/>
    SELECT ?label
    from <http://purl.obolibrary.org/obo/merged/GO>
    WHERE
    {
       ?annotation rdfs:label ?label .
       VALUES ?annotation $terms .
    }
    """
end

import inspect
from pysb import *
import pysb.core
from pysb.core import ComponentSet, as_reaction_pattern, as_complex_pattern, MonomerPattern, ComplexPattern
import numbers
import functools
import itertools

# Internal helper functions (Copied from PySB)
# =========================

def _complex_pattern_label(cp):
    """Return a string label for a ComplexPattern."""
    mp_labels = [_monomer_pattern_label(mp) for mp in cp.monomer_patterns]
    return ''.join(mp_labels)

def _monomer_pattern_label(mp):
    """Return a string label for a MonomerPattern."""
    site_values = [str(x) for x in mp.site_conditions.values()
                            if x is not None
                            and not isinstance(x, list)
                            and not isinstance(x, tuple)
                            and not isinstance(x, numbers.Real)]
    return mp.monomer.name + ''.join(site_values)

def _rule_name_generic(rule_expression):
    """Return a generic string label for a RuleExpression."""
    # Get ReactionPatterns
    react_p = rule_expression.reactant_pattern
    prod_p = rule_expression.product_pattern
    # Build the label components
    lhs_label = [_complex_pattern_label(cp) for cp in react_p.complex_patterns]
    lhs_label = '_'.join(lhs_label)
    rhs_label = [_complex_pattern_label(cp) for cp in prod_p.complex_patterns]
    rhs_label = '_'.join(rhs_label)
    return '%s_to_%s' % (lhs_label, rhs_label)

def _macro_rule(rule_prefix, rule_expression, klist, ksuffixes,
                name_func=_rule_name_generic):
    """
    A helper function for writing macros that generates a single rule.

    Parameters
    ----------
    rule_prefix : string
        The prefix that is prepended to the (automatically generated) name for
        the rule.
    rule_expression : RuleExpression
        An expression specifying the form of the rule; gets passed directly
        to the Rule constructor.
    klist : list of Parameters or Expressions, or list of numbers
        If the rule is unidirectional, the list must contain one element
        (either a Parameter/Expression or number); if the rule is reversible,
        it must contain two elements. If the rule is reversible, the first
        element in the list is taken to be the forward rate, and the second
        element is taken as the reverse rate. 
    ksuffixes : list of strings
        If klist contains numbers rather than Parameters or Expressions, the
        strings in ksuffixes are used to automatically generate the necessary
        Parameter objects. The suffixes are appended to the rule name to
        generate the associated parameter name. ksuffixes must contain one
        element if the rule is unidirectional, two if it is reversible.
    name_func : function, optional
        A function which takes a RuleExpression and returns a string label for
        it, to be called as part of the automatic rule name generation. If not
        provided, a built-in default naming function will be used.

    Returns
    -------
    components : ComponentSet
        The generated components. Contains the generated Rule and up to two
        generated Parameter objects (if klist was given as numbers).

    Notes
    -----
    The default naming scheme (if `name_func` is not passed) follows the form::

        '%s_%s_to_%s' % (rule_prefix, lhs_label, rhs_label)

    where lhs_label and rhs_label are each concatenations of the Monomer names
    and specified sites in the ComplexPatterns on each side of the
    RuleExpression. The actual implementation is in the function
    _rule_name_generic, which in turn calls _complex_pattern_label and
    _monomer_pattern_label. For some specialized reactions it may be helpful to
    devise a custom naming scheme rather than rely on this default.

    Examples
    --------
    Using distinct Monomers for substrate and product::

        >>> from pysb import *
        >>> from pysb.macros import _macro_rule
        >>> 
        >>> Model() # doctest:+ELLIPSIS
        <Model '_interactive_' (monomers: 0, rules: 0, parameters: 0, expressions: 0, compartments: 0) at ...>
        >>> Monomer('A', ['s'])
        Monomer('A', ['s'])
        >>> Monomer('B', ['s'])
        Monomer('B', ['s'])
        >>> 
        >>> _macro_rule('bind', A(s=None) + B(s=None) <> A(s=1) % B(s=1),
        ... [1e6, 1e-1], ['kf', 'kr']) # doctest:+NORMALIZE_WHITESPACE
        ComponentSet([
         Rule('bind_A_B_to_AB', A(s=None) + B(s=None) <> A(s=1) % B(s=1),
             bind_A_B_to_AB_kf, bind_A_B_to_AB_kr),
         Parameter('bind_A_B_to_AB_kf', 1000000.0),
         Parameter('bind_A_B_to_AB_kr', 0.1),
         ])

    """

    r_name = '%s_%s' % (rule_prefix, name_func(rule_expression))

    # If rule is unidirectional, make sure we only have one parameter
    if (not rule_expression.is_reversible):
        if len(klist) != 1 or len(ksuffixes) != 1:
            raise ValueError("A unidirectional rule must have one parameter.")
    # If rule is bidirectional, make sure we have two parameters
    else:
        if len(klist) != 2 or len(ksuffixes) != 2:
            raise ValueError("A bidirectional rule must have two parameters.")

    if all(isinstance(x, (Parameter, Expression)) for x in klist):
        k1 = klist[0]
        if rule_expression.is_reversible:
            k2 = klist[1]
        params_created = ComponentSet()
    # if klist is numbers, generate the Parameters
    elif all(isinstance(x, numbers.Real) for x in klist):
        k1 = Parameter('%s_%s' % (r_name, ksuffixes[0]), klist[0])
        params_created = ComponentSet([k1]) 
        if rule_expression.is_reversible:
            k2 = Parameter('%s_%s' % (r_name, ksuffixes[1]),
                           klist[1])
            params_created.add(k2)
    else:
        raise ValueError("klist must contain Parameters, Expressions, or numbers.")

    if rule_expression.is_reversible:
        r = Rule(r_name, rule_expression, k1, k2)
    else:
        r = Rule(r_name, rule_expression, k1)

    # Build a set of components that were created
    return ComponentSet([r]) | params_created

def _verify_sites(m, *site_list):
    """
    Checks that the monomer m contains all of the sites in site_list.

    Parameters
    ----------
    m : Monomer or MonomerPattern
        The monomer to check.
    site1, site2, ... : string
        One or more site names to check on m

    Returns
    -------
    True if m contains all sites; raises a ValueError otherwise.

    Raises
    ------
    ValueError
        If any of the sites are not found.

    """

    if isinstance(m, ComplexPattern):
        return _verify_sites_complex(m, *site_list)
    else:
        for site in site_list:
            if site not in m().monomer.sites:
                raise ValueError("Monomer '%s' must contain the site '%s'" %
                                (m().monomer.name, site))
        return True

def _verify_sites_complex(c, *site_list):
  
    """
    Checks that the complex c contains all of the sites in site_list.

    Parameters
    ----------
    c : ComplexPattern
        The complex to check.
    site1, site2, ... : string
        One or more site names to check on c

    Returns
    -------
    If all sites are found within the complex, a dictionary of monomers and the sites within site_list they contain.  Raises a ValueError if one or more sites not in the complex.

    Raises
    ------
    ValueError
         If any of the sites are not found within the complex.

    """

    allsitesdict = {}
    for mon in c.monomer_patterns:
        allsitesdict[mon] = mon.monomer.sites
    for site in site_list:
        specsitesdict = {}
        for monomer, li in allsitesdict.items():
            for s in li:
                if site in li:
                    specsitesdict[monomer] = site
        if len(specsitesdict) == 0:
            raise ValueError("Site '%s' not found in complex '%s'" % (site, c))
    return specsitesdict

def bind_name_func(rule_expression):
    react_cps = rule_expression.reactant_pattern.complex_patterns
    return '_'.join(_complex_pattern_label(cp) for cp in react_cps)

# Boolean to Rules
# ================

# def macro_monomer(monomer_prefix, site_list, m_0, init_states, site_states=None): # possible helper function
#     """
#     A helper function for writing macros that generates a single monomer.
# 
#     Parameters
#     ----------
#     monomer_prefix : string
#         The prefix that is prepended to the (automatically generated) name for
#         the monomer.
#     site_list : list of monomer sites
#     m_0 : numerical value
#         Initial value of the monomer
#     init_states: list of initial site states
#     site_states: dict of possible site states
# 
#     Returns
#     -------
#     components : ComponentSet
#         The generated components. Contains the generated Monomer, Initial condition parameter
#         and observables for each possible state (if bound/unbound)
# 
# 
#     """
#     m_name = '%s' % (monomer_prefix)
#     m = Monomer(m_name, site_list, site_states)
#     dict = {site_list[0]:None, site_list[1]:None}
#     m_Pattern = MonomerPattern(monomer=m, site_conditions=dict, compartment=None)
#     p = Parameter('%s_%s' % (m_name, '0'), m_0)
#     Initial(m_Pattern, p)
#     
#     components = [m,p]
#     if site_states == None:
#         
#         for each in itertools.product('10', repeat = len(init_states)):
#             combo = ''.join(each)
#     #         print combo
#             statelist = [j for j in range(len(each))]
#             for j in range(len(each)):
#                 if each[j] == '0':
#                     statelist[j] = None
#             if len(init_states) == 2:
#                 obs = {site_list[0]:statelist[0], site_list[1]:statelist[1]}
#                 o_Pattern = MonomerPattern(monomer=m, site_conditions=obs, compartment=None)
#                 components.append(Observable('%s_%s_%s' % (m_name, combo, 'obs'), o_Pattern))
# 
#     return ComponentSet(components)
    
# Transcription

def transcription_switch(gene, g_sites, TF_list, TF_sites, klist, product, tlist, deg_r, switch_state):
    
    """
    Generates a set of reactions describing transcriptional regulation. These reactions 
    account for both the mass action aspects of trascription factor binding and the 
    combinatorial/switch-like apects of transcription
    
    Parameters
    ----------
    gene : Monomer
    TF_list: list of Monomers
        List of transcription factors, both activators and inhibitors that may bind the gene
    g_sites, TF_sites : string
        The names of the binding sites on the gene and TFs respectively.
    klist : list of 2*|TF_list| Parameters
        The forward and reverse rate constants for each TF.
    tlist : list of 2 Parameters
        list of 2 rate constants, an active rate and a leakage rate
    product : Monomer
    switch_state : string - binary
        A Boolean state describing the 'active' configuration of bound transcription
        factors
    
    Returns
    -------
    components : ComponentSet
        The generated components.
    """
       
    components = ComponentSet()
    
    for i in range(len(TF_list)):
        if isinstance(TF_list[i], list):
            for j in range(len(TF_list[i])):
                components |= _macro_rule('bind'+str(i),
                                    gene({g_sites[i]: None}) + TF_list[i][j]({TF_sites[i]: None}) <>
                                    gene({g_sites[i]: 1}) % TF_list[i][j]({TF_sites[i]: 1}),
                                    klist[i*2:i*2+2], ['kf', 'kr'], name_func=bind_name_func)
        else:
            components |= _macro_rule('bind'+str(i),
                    gene({g_sites[i]: None}) + TF_list[i]({TF_sites[i]: None}) <>
                    gene({g_sites[i]: 1}) % TF_list[i]({TF_sites[i]: 1}),
                    klist[i*2:i*2+2], ['kf', 'kr'], name_func=bind_name_func)              

    for each in itertools.product('10', repeat = len(g_sites)):
        combo = ''.join(each)
        statelist = []
        for j in range(len(each)):
            if each[j] == '1':
                statelist.append(j)
            else:
                statelist.append(None)

        site_states = {}
        for j in range(len(g_sites)):
            site_states[g_sites[j]] = statelist[j]
#         print site_states

        if combo == switch_state:
            components |= _macro_rule('transcription'+str(combo),
                                gene(site_states) >>
                                gene(site_states) + product(),
                                [tlist[0]], ['kt'], name_func=bind_name_func)
        else:
            components |= _macro_rule('transcription'+str(combo),
                                gene(site_states) >>
                                gene(site_states) + product(),
                                [tlist[1]], ['kl'], name_func=bind_name_func)
                
    components |= _macro_rule('degradation', product() >> None, [deg_r], ['kd'], name_func=bind_name_func)
   
    return components

def transcription(gene, g_sites, TF_list, TF_sites, klist, product, tlist, deg_r):
    """
    Generates a set of reactions describing the transcriptional regulation. These reactions 
    account for both the mass action aspects of trascription factor binding and the 
    combinatorial/switch-like apects of transcription
    
    Parameters
    ----------
    gene : Monomer
    TF_list: list of Monomers
        List of transcription factors, both activators and inhibitors that may bind the gene
    g_sites, TF_sites : string
        The names of the binding sites on the gene and TFs respectively.
    klist : list of 2*|TF_list| Parameters
        The forward and reverse rate constants for each TF.
    tlist : list of 2^|TF_list| Parameters
        list of rate constants, one for each combination of TFs
    product : Monomer
    
    Returns
    -------
    components : ComponentSet
        The generated components.
    """
       
    components = ComponentSet()
    
    for i in range(len(TF_list)):
        if isinstance(TF_list[i], list):
            for j in range(len(TF_list[i])):
                components |= _macro_rule('bind'+str(i),
                                    gene({g_sites[i]: None}) + TF_list[i][j]({TF_sites[i]: None}) <>
                                    gene({g_sites[i]: 1}) % TF_list[i][j]({TF_sites[i]: 1}),
                                    klist[i*2:i*2+2], ['kf', 'kr'], name_func=bind_name_func)
        else:
            components |= _macro_rule('bind'+str(i),
                    gene({g_sites[i]: None}) + TF_list[i]({TF_sites[i]: None}) <>
                    gene({g_sites[i]: 1}) % TF_list[i]({TF_sites[i]: 1}),
                    klist[i*2:i*2+2], ['kf', 'kr'], name_func=bind_name_func)                    
    
    for i,each in enumerate(itertools.product('01', repeat = len(g_sites))):
        combo = ''.join(each)
        statelist = []
        for j in range(len(each)):
            if each[j] == '1':
                statelist.append(j)
            else:
                statelist.append(None)
                
        site_states = {}
        for j in range(len(g_sites)):
            site_states[g_sites[j]] = statelist[j]
              
        components |= _macro_rule('transcription'+str(combo),
                            gene(site_states) >>
                            gene(site_states) + product(),
                            [tlist[i]], [combo], name_func=bind_name_func)

    components |= _macro_rule('degradation', product() >> None, [deg_r], ['kd'], name_func=bind_name_func)
    
    return components

def translation(mRNA, protein, klist):
    """
    Generates a set of reactions describing the translation of mRNA to protein.
    
    Parameters
    ----------
    mRNA : Monomer
    k_list : list of 2 Parameters
        List of two rate constants for translation and degradation of protein.
    product : Monomer
    Returns
    -------
    components : ComponentSet
        The generated components.
    """

    site_states = {}
    for j in range(len(protein.sites)):
        site_states[protein.sites[j]] = None

    components = _macro_rule('translation',
                        mRNA() >> mRNA() + protein(site_states), 
                        [klist[0]], ['kt'], name_func=bind_name_func)
     
    components |= _macro_rule('degradation', 
                              protein(site_states) >> None, [klist[1]], ['kd'], name_func=bind_name_func)
  
    return components

def translation_persistant(mRNA, protein, klist):
    """
    Generates a reaction describing the translation of mRNA to protein.
    NO degradation reaction is provided
    
    Parameters
    ----------
    mRNA : Monomer
    k_list : list of 2 Parameters
        List of two rate constants for translation and degradation of protein.
    product : Monomer
    Returns
    -------
    components : ComponentSet
        The generated components.
    """

    site_states = {}
    for j in range(len(protein.sites)):
        site_states[protein.sites[j]] = None

    components = _macro_rule('translation',
                        mRNA() >> mRNA() + protein(site_states), 
                        [klist[0]], ['kt'], name_func=bind_name_func)
  
    return components

def binary_transformation(enzyme, e_site, substrate, s_site, product_1, product_2, k_list):
    """
    Generates a set of reactions describing a catalytic reaction that results in the consumption 
    of substrate and subsequent production of one of two products depending on the presence of enzyme. 
    (could be expanded to accept additional products and/or enzymes). Additional sites are assumed 
    to be unbound.
    
    Parameters
    ----------
    enzyme, substrate : Monomers
    substrate : Monomer
    e_site : List of strings
        binding and activation sites on the enzyme
    s_site : String
        binding site on the substrate
    product_1, product_2 : Monomers
    
    k_list : list of 6 Parameters
        parameters represent, in the following order,
        - the forward and reverse binding rates for enzyme and substrate
        - the catalytic rate constant for product_1 (with enzyme present)
        - the catalytic rate constant for product_2 (without enzyme present)
        - the degradation rates for product_1 and product_2
        
    Returns
    -------
    components : ComponentSet
        The generated components.
    """
    
    enzyme_sites_pre = {}
    for each in enzyme.sites:
        enzyme_sites_pre[each] = None

    enzyme_sites_bound = {}
    for each in enzyme.sites:
        if each == e_site:
            enzyme_sites_bound[each] = 1
        else:
            enzyme_sites_bound[each] = None
            
    components = ComponentSet()
            
    components |= _macro_rule('bind',
            enzyme(enzyme_sites_pre) + substrate({s_site: None}) <>
            enzyme(enzyme_sites_bound) % substrate({s_site: 1}),
            [k_list[0],k_list[1]], ['kf', 'kr'], name_func=bind_name_func)  

    pro_sites_1 = {}
    for each in product_1.sites:
        pro_sites_1[each] = None  

    components |= _macro_rule('catalyze_1',
            enzyme(enzyme_sites_bound) % substrate({s_site: 1}) >>
            enzyme(enzyme_sites_pre) + product_1(pro_sites_1),
            [k_list[2]], ['kc1'], name_func=bind_name_func)

    pro_sites_2 = {}
    for each in product_2.sites:
        pro_sites_2[each] = None 

    components |= _macro_rule('catalyze_2',
            substrate({s_site: None}) >> product_2(pro_sites_2),
            [k_list[3]], ['kc2'], name_func=bind_name_func)
    
    components |= _macro_rule('degradation_1', product_1() >> None, [k_list[4]], ['kd1'], name_func=bind_name_func)
    components |= _macro_rule('degradation_2', product_2() >> None, [k_list[5]], ['kd2'], name_func=bind_name_func)
    
    return components
    
def dimerization(s1, site_1, s2, site_2, product, k_list):
    """
    This is essentially identical to a binding reaction. The differences are the
    possibility of multiple binding partners and the relabeling as a new species 
    to be consistent with the logical model.
    
    Currently assumes Monomers are completely unbound
    
    Parameters
    ----------
    s1 : Monomer
    s2 : List of identicle Monomers
        Monomers participating in the binding reaction.
    site_1, site_2 : string, list of strings
        The names of the sites on s1 and s2 used for binding.
    product : 
    k_list : list of 3 Parameters
        Forward and reverse rate constants (in that order) and 
        a degradation rate for the product.

    Returns
    -------
    components : ComponentSet
        The generated components. 
    """
    
    pro_sites = {}
    for each in product.sites:
        pro_sites[each] = None
    
    s1_sites = {}
    for each in s1.sites:
        if each == site_1:
            s1_sites[each] = None
        else:
            s1_sites[each] = None
    
    s2_sites = {}
    for each in s2[0].sites:
        if each == site_2:
            s2_sites[each] = None
        else:
            s2_sites[each] = None
    
    components = ComponentSet()
    
    for each in s2:
        components |= _macro_rule('dimerize',
                           s1(s1_sites) + each(s2_sites) <> product(pro_sites),
                           k_list[0:2], ['kf', 'kr'], name_func=bind_name_func)
    
    components |= _macro_rule('degradation', product() >> None, [k_list[2]], ['kd'], name_func=bind_name_func)

def proteolysis(enzyme, e_site, target, t_site, k_list):
    """
    This represents degradation of one species by another.
    
    Parameters
    ----------
    enzyme, target : Monomers
    e_site, t_site : strings
        The names of the sites on s1 and s2 used for binding.
    product : 
    k_list : list of 3 Parameters
        Forward and reverse rate constants (in that order) and 
        a degradation rate.

    Returns
    -------
    components : ComponentSet
        The generated components. 
    """
    
    components = ComponentSet()
    
    components |= _macro_rule('bind',
                       enzyme({e_site:None}) + target({t_site:None}) <> enzyme({e_site:1}) % target({t_site:1}),
                       k_list[0:2], ['kf', 'kr'], name_func=bind_name_func)
    
    components |= _macro_rule('proteolysis', enzyme({e_site:1}) % target({t_site:1}) >> enzyme({e_site:None}), [k_list[2]], ['kd'], name_func=bind_name_func)



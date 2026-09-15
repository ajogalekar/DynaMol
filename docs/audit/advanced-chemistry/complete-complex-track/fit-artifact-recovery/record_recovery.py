"""Record bounded primary-source retrieval; no chemistry parameters assigned."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
CACHE=Path('/Users/ashujo/.cache/dynamol-research/complete-complex-track/fit-artifact-recovery-v1')
TRACK=HERE.parent
QUERIES=[
    '"Kenneth J. Miller" "Additivity" polarizability pdf rpi',
    'site:q4md-forcefieldtools.org "dimethyl" "phosphate"',
    'site:theory.rutgers.edu "Panteva" parameters download DMP',
    'site:github.com "c6" "lj_1264_pol"',
    '"Additivity Methods in Molecular Polarizability" filetype:pdf -researchgate -electronicsandbooks',
    '"Miller" "8533" "CTE" polarizability',
    '"DMP" "Dupradeau" mol2',
    '"dimethyl phosphate" "R.E.DD.B."',
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    prior_paths=[TRACK/'c6-gdp-transfer/REPORT.md',TRACK/'c6-gdp-transfer/evidence.json',
                 TRACK/'PANTEVA-APPLICABILITY.md',TRACK/'panteva-applicability.json']
    hashes={str(p):sha(p) for p in prior_paths}
    source=json.loads((TRACK/'c6-gdp-transfer/evidence.json').read_text())
    retained_scientific={key:value for key,value in source['source_hashes'].items()
                         if key in {'native_topology','mapping','assembly','gdp_scoping','gdp_source_prep','metal_panel','covalent_panel'}}
    for row in retained_scientific.values():
        if sha(row['path'])!=row['sha256']:
            raise ValueError('Retained molecular source changed')
    initial=json.loads((CACHE/'retrieval-initial.json').read_text())
    final=json.loads((CACHE/'retrieval-final.json').read_text())
    projects=json.loads((CACHE/'retrieval-projects.json').read_text())
    download=json.loads((CACHE/'retrieval-download.json').read_text())
    metadata=json.loads((CACHE/'miller-oa-metadata.json').read_text())
    db=(CACHE/'red-db-institutional.html').read_text()
    form=(CACHE/'red-db-download.html').read_text()
    pub=(CACHE/'york-publication.html').read_text()
    if not all(text in form for text in ['method="POST" action="download.php"','value = "bymolname"','name="enter"','name="ok" value="OK"']):
        raise ValueError('Observed public database form differs from recorded next step')
    source_files={str(p):{'sha256':sha(p),'bytes':p.stat().st_size} for p in CACHE.iterdir() if p.is_file()}
    result={
        'schema_version':1,'stage':'bounded_fit_artifact_recovery','recorded_utc':datetime.now(timezone.utc).isoformat(),
        'limits':{'targeted_search_queries_maximum':8,'targeted_search_queries_used':len(QUERIES),'elapsed_limit_minutes':10,
                  'scope':'Targeted search-query count excludes ordinary retrieval of discovered source pages and public archive metadata. No search form was submitted after reaching the eight-query limit.'},
        'queries':QUERIES,
        'prior_evidence_sha256':hashes,'retained_scientific_sources':retained_scientific,
        'source_files':source_files,'retrievals':initial+final+[projects,download],
        'web_tool_errors':[{'url':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2918240/fullTextXML','result':'Web reader reported nonretryable URL-open failure; independent ordinary HTTPS retrieval returned404.'},
                           {'url':'https://q4md-forcefieldtools.org/REDDB/index.php','result':'Web reader timed out; ordinary HTTPS retrieval failed certificate hostname verification. Verification was not disabled.'}],
        'new_primary_route':{
            'database':'RESP ESP charge DDataBase (R.E.DD.B.)',
            'institutional_base':'https://upjv.q4md-forcefieldtools.org/REDDB/',
            'host_discovery':'Indexed author-hosted q4md document used the upjv subdomain; retrieved database links and author branding independently confirm that host.',
            'TLS_verification_enabled':True,'index_accessible_without_authentication':True,
            'list_projects_accessible':True,'public_search_download_form_accessible':True,
            'index_sha256':sha(CACHE/'red-db-institutional.html'),
            'first_listing_scope':'First default page only:10 fragment projects. Not a complete project search and not evidence that DMP is absent.',
            'next_observed_request':{'method':'POST','url':'https://upjv.q4md-forcefieldtools.org/REDDB/download.php',
                                     'fields':{'choice':'bymolname','enter':'dimethyl phosphate','ok':'OK'},'executed':False},
            'dataset_license_established':False,
            'license_note':'Footer labels projects free and links to GNU philosophy; no project-specific redistribution license or exact Panteva fit provenance was established.'},
        'Miller':{'doi':'10.1021/ja00179a044','original_table_recovered':False,
                  'metadata_discovery_only':{'source':'OpenAlex','source_is_not_primary_scientific_evidence':True,
                                            'open_access':metadata.get('open_access'),'best_oa_location':metadata.get('best_oa_location'),
                                            'locations':metadata.get('locations')},
                  'institutional_author_archive_found':False,
                  'search_result_mirror_limit':'Search results surfaced an unrelated scanned-journal mirror; no institution/author provenance or table extraction was established, and it was not used for a c6 mapping.',
                  'documented_c6_crosswalk_recovered':False,'alpha_assigned':None},
        'DMP':{'exact_charge_or_topology_artifact_recovered':False,'database_candidate_project_identified':False,
               'new_accessible_primary_database_route_identified':True,
               'exact_link_to_Panteva_fit_established':False,
               'numerical_DMP_GDP_charge_comparison_performed':False,
               'reason_to_stop_this_pass':'Eight targeted search queries completed. A newly accessible public database form is recorded for the next bounded lookup, rather than exceeding this pass\'s query cap.'},
        'Panteva_author_page':{'url':'https://theory.rutgers.edu/publications_page.php?publication_id=25',
                               'sha256':sha(CACHE/'york-publication.html'),
                               'linked_article_or_parameter_like_urls':sorted(set(re.findall(r'href=[\"\']([^\"\']+(?:pdf|mol2|prmtop|zip|tgz|tar\.gz))',pub))),
                               'finding':'Visible full-article and paper-PDF links; no named DMP topology/charge bundle was found on this publication page.'},
        'next_bounded_step':'Submit the observed public R.E.DD.B. molecule-name form once for dimethyl phosphate. Inspect returned project metadata, charge/model version, files and rights; do not equate a matching molecule/project with the exact Panteva2015 fit without a provenance link.',
        'not_a_stop_of_overall_advanced_track':True,
        'new_md_or_qm':False,'parameters_assigned_or_changed':False,'topology_changed':False,
        'app_or_shared_checkpoint_changed':False,'paid_services_used':False,'authentication_bypassed':False,
        'messages_to_authors_sent':False,'app_ready':False,'physical_model_validated':False}
    if hashes!={str(p):sha(p) for p in prior_paths}:
        raise ValueError('Prior report changed while recording recovery')
    result['all_retained_scientific_sources_unchanged']=True
    (HERE/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'queries':len(QUERIES),'new_primary_database_route':True,'exact_fit_artifact_recovered':False,
                      'c6_crosswalk_recovered':False,'sources_unchanged':True,'output':str(HERE/'evidence.json')},indent=2))


if __name__=='__main__':
    main()

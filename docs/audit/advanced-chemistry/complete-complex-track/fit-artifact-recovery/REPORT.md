# Bounded fit-artifact recovery — 14 September 2026

**A working primary database route was recovered, but no exact DMP fit artifact or `c6` crosswalk was recovered in this pass.** Eight targeted search queries were completed. No parameters, topology, or application behavior changed.

The [evidence JSON](evidence.json) retains the query list, URLs, source hashes, failures, and observed database form. [record_recovery.py](record_recovery.py) verifies the existing scientific inputs and records the result without molecular calculations.

The useful new route is the author-maintained [R.E.DD.B. institutional host](https://upjv.q4md-forcefieldtools.org/REDDB/). Its index, first project-list page, and public download/search form load successfully with certificate verification enabled and no login. The bare `q4md-forcefieldtools.org` hostname fails certificate hostname verification; verification was not disabled. The working subdomain was discovered in an indexed author-hosted document and confirmed by the retrieved database itself.

The [public download page](https://upjv.q4md-forcefieldtools.org/REDDB/download.php) provides an observed molecule-name search form. Its next request is recorded as a POST with `choice=bymolname`, `enter=dimethyl phosphate`, and `ok=OK`. It was not submitted after reaching this pass's eight-query cap. Only the first default project-list page was inspected; it does not establish whether DMP is present or absent. The footer's description of projects as “free” is not treated as a verified project-specific redistribution license.

Other genuinely new routes did not supply the missing artifact:

- Europe PMC's public full-text XML route for the Dupradeau paper and the attempted NCBI open-access metadata route returned 404. Previously challenged publisher/PMC pages were not retried.
- The [Rutgers author publication page](https://theory.rutgers.edu/publications_page.php?publication_id=25) supplies the article and paper-PDF links, but no named DMP charge/topology bundle was found there.
- OpenAlex discovery metadata lists Miller's 1990 article as closed, with no repository full-text location. That metadata is used solely for locating sources, not as scientific evidence. An unrelated scanned-journal mirror appeared in search results, but no author/institution provenance or table extraction was established; it was not used to infer a `c6` value.

All prior evidence and the retained topology, source mapping, GDP model, and fixed case panels remain unchanged. No exact project, partial-charge set, or publication link has yet been established for the DMP model used in the Panteva fit. Attribution to Dupradeau alone still does not identify that artifact. The defined baseline GDP O2/O3/NB coverage from the preceding audit remains distinct from transferring fitted corrections.

Research downloads and retrieval logs are retained in:

`/Users/ashujo/.cache/dynamol-research/complete-complex-track/fit-artifact-recovery-v1/`

The next narrow lookup is the observed public R.E.DD.B. molecule-name form. Any returned DMP project needs its authorship, date, molecular charge, file version, rights, and explicit relationship to the Panteva calculation checked before it can be called the original fit artifact. This ends the bounded retrieval pass, not the separate advanced-chemistry track.

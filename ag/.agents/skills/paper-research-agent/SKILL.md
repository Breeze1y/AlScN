---
name: paper-research-agent
description: Search primary research papers for this Ag/AlN/AlScN/AlN/Pt memristor project, verify sources, and turn findings into evidence-based Chinese summaries with concrete implications for the current model and experiments. Trigger when asked to find, update, compare, or summarize literature for this project.
---

# Project paper research agent

Use this workflow when the user asks for literature search, a literature review, paper comparisons, or project-relevant summaries in this repository.

## 1. Read the current project before searching

Read `README.md`, `docs/毕业设计完整执行方案-20261104.md`, and the most relevant current report or model documentation. Inspect `papers/` to see what has already been curated. Treat the files and dates in the repository as the source of truth for current implementation status; do not repeat old project summaries as if they were current.

Build a short working brief that distinguishes:

- the device stack and experimental protocols;
- implemented code and validated evidence;
- exploratory assumptions and unvalidated parameters;
- explicit limitations and pending work.

This repository studies an Ag/AlN/AlScN/AlN/Pt resistive memory device using effective KMC and electrical models. Verify the exact current scope from repository files on every run.

## 2. Search in focused passes

Translate the project question into a few searchable themes. For this project, useful terms include `CBRAM`, `ECM`, `Ag filament`, `AlN`, `AlScN`, `kinetic Monte Carlo`, `filament morphology`, `RESET rupture`, `retention`, `polarization`, `interface barrier`, `microvoid`, and `residual stress`.

Search first for the exact material/stack and then for adjacent material systems or methods. Use multiple concise English queries and, where useful, Chinese equivalents. Search recent work as well as foundational papers. Follow citation trails only after identifying the original paper.

Prioritize original research papers and primary sources such as journal publishers, conference publishers, PubMed records, and author or institutional repositories. Reviews can help find references, but do not treat them as primary evidence. Open the source and check its abstract, method, results, and limitations where available. Verify title, authors, year, venue, and DOI against a publisher, DOI, or institutional record.

If the full text is paywalled or unavailable, say which parts were actually accessible. Never imply that a paper was fully read based only on a search snippet. Use DOI links in the report and cite every substantive claim close to where it appears.

## 3. Rate relevance and evidence

Classify each selected paper as:

1. **Direct device match** — same or nearly identical stack and switching mechanism.
2. **Material-adjacent evidence** — e.g. Ag migration in AlN or AlScN polarization, but a different device or process.
3. **Method-adjacent evidence** — KMC, tomography, electrical extraction, or reliability methods from another material system.

For each paper, record the device/material, method, observed result, evidence source, relevance to this project, and a concrete transfer limit. Separate:

- **Paper finding**: what the authors measured or simulated.
- **Project implication**: a proposed use in this repository, explicitly labeled as an inference when it is not directly demonstrated by the source.

Never transfer a migration barrier, reaction barrier, permittivity, polarization, current density, or time scale from another material stack as a measured Ag/AlN/AlScN value. Do not describe a proposed mechanism as established for this device unless the paper or project data directly support it.

## 4. Connect literature to the implementation

Map findings to existing components, model assumptions, and planned work. Check whether the project already has the proposed feature before recommending it. In particular, distinguish a fixed effective polarization field from dynamic ferroelectric switching; distinguish first-bridge forming from a complete SET/RESET/retention cycle; distinguish a two-dimensional local window from the full device geometry; and distinguish geometric connectivity from electrical HRS/LRS labels.

Recommendations should name the observable or comparison that would make the literature useful, such as bridge probability, forming-voltage distribution, filament width/roughness, minimum cross-section, reset rupture class, read resistance, retention drift, or cycle-to-cycle spread. Prefer small controlled comparisons and shared protocols over adding unvalidated physical parameters.

Do not modify model code, change calibrated parameters, download or redistribute paywalled PDFs, or claim the project has reproduced a paper unless the user explicitly asks for that work.

## 5. Write a project-ready report

For a substantial literature task, save a dated Markdown report under `reports/literature/`. Use `references/report-template.md` as the default structure. The report should include:

- a concise synthesis for the current research question;
- the project state and assumptions checked, with repository file links;
- a table of selected papers with DOI/publisher links and relevance classes;
- what each source supports and what it does not support;
- prioritized next steps tied to model observables or existing datasets;
- the search date, query themes, source-access limits, and open questions.

For a narrow request, answer in chat and only create a file if it would materially help the project. Keep quotations short; paraphrase and cite. Avoid duplicating existing paper summaries: refer to local `papers/` files and clearly identify newly found literature.

## 6. Handoff

Finish with the report path, the most useful findings, and any material source-access limitation. If the search did not find direct evidence, say so plainly and summarize the closest evidence without calling it a direct match.

// Offline unit tests for the isolate-validator helpers added in v1.21.13:
// lineage/species parsing, the coverage -> IUPAC degeneracy recommender, and the
// summary / lineage / coverage renderers. Boots the real index.html <script> in a
// stubbed context (same pattern as ui_workbench.js) and calls the functions directly.
const fs = require("fs"), vm = require("vm");
const script = fs.readFileSync("static/index.html", "utf8").match(/<script>([\s\S]*)<\/script>/)[1];
const cache = {};
function el(id){ return cache[id] || (cache[id] = {_id:id, value:"", checked:false, innerHTML:"", textContent:"",
  style:{}, classList:{add(){},remove(){},contains:()=>false}, appendChild(){}, addEventListener(){}, click(){},
  querySelectorAll:()=>[], querySelector:()=>({value:""}) }); }
const store={};
const ctx={ document:{getElementById:el, createElement:()=>el("c"+Math.random()),
    querySelectorAll:(s)=>(s==="nav a"||s==="section")?Array.from({length:16},(_,i)=>el("n"+i)):[]},
  window:{}, console, Blob:function(){}, URL:{createObjectURL:()=>"x",revokeObjectURL(){}},
  fetch:async()=>({ok:true,status:200,json:async()=>({})}),
  localStorage:{getItem:k=>store[k]||null,setItem:(k,v)=>{store[k]=v;}}, alert:()=>{}, confirm:()=>true,
  setTimeout, JSON, Math, Object, Array, Date, parseInt, parseFloat, isFinite, Number, String, Boolean, RegExp };
vm.createContext(ctx); vm.runInContext(script, ctx); ctx.api=async()=>({});

let ok=0, bad=[];
const T=(label,fn)=>{ try{ fn(); ok++; console.log("  \u2713 "+label); }catch(e){ bad.push(label); console.log("  \u2717 "+label+" -> "+e.message); } };
const eq=(a,b,m)=>{ if(a!==b) throw new Error((m||"")+" expected "+JSON.stringify(b)+" got "+JSON.stringify(a)); };
const has=(s,sub,m)=>{ if((s||"").indexOf(sub)<0) throw new Error((m||"")+" missing "+JSON.stringify(sub)); };

// --- isoLineage: species + lineage token parsing on real-world titles ---
T("isoLineage: named species -> 'P. relictum'", ()=>{
  eq(ctx.isoLineage("Plasmodium relictum strain SGS1 genome assembly, chromosome: 1").species, "P. relictum"); });
T("isoLineage: 'sp.' record -> genus + lineage token", ()=>{
  const l=ctx.isoLineage("UNVERIFIED: Plasmodium sp. isolate TX21010_r1 cytochrome b (cytb) gene, partial sequence");
  eq(l.species,"Plasmodium","species"); eq(l.lineage,"TX21010_r1","lineage"); });
T("isoLineage: haplotype token", ()=>{
  eq(ctx.isoLineage("Plasmodium sp. haplotype MALCIN01 cytochrome b").lineage, "MALCIN01"); });
T("isoLineage: different genus", ()=>{
  eq(ctx.isoLineage("Haemoproteus majoris isolate WW2 cytochrome b").species, "H. majoris"); });

// --- IUPAC helpers ---
T("iuCode / iuMatch", ()=>{
  eq(ctx.iuCode(["A","G"]),"R"); eq(ctx.iuCode(["C","T"]),"Y"); eq(ctx.iuCode(["A","C","G","T"]),"N");
  if(!ctx.iuMatch("R","A")) throw new Error("R should match A");
  if(ctx.iuMatch("R","C")) throw new Error("R should NOT match C");
  if(!ctx.iuMatch("W","T")) throw new Error("W should match T"); });
T("stripMods: IDT order string -> bare", ()=>{
  eq(ctx.stripMods("/56-FAM/CTTA+CA+A+GATAT+CC+ACCACA/3IABkFQ/"), "CTTACAAGATATCCACCACA"); });

// --- isoOligoCoverage: per-base conservation + degeneracy suggestion ---
const F="ACGTACGTAC";                       // idx5 = 'C'
const mk=(w)=>({role:"target", f_win:w, r_win:null, p_win:null});
T("coverage: minor allele in >=2 isolates -> suggest IUPAC", ()=>{
  const rows=[mk(F),mk(F),mk(F),mk("ACGTATGTAC"),mk("ACGTATGTAC")];  // pos6 C->T in 2/5
  const c=ctx.isoOligoCoverage(rows, F, "f_win");
  eq(c.n,5,"n"); eq(c.suggested,"ACGTAYGTAC","suggested"); eq(c.nchg,1,"nchg");
  eq(c.variants.length,1,"variants"); eq(c.variants[0].k,6,"pos"); eq(c.variants[0].pct,60,"pct"); });
T("coverage: singleton allele does NOT trigger a suggestion", ()=>{
  const rows=[mk(F),mk(F),mk(F),mk(F),mk("ACGTATGTAC")];             // T in only 1/5
  const c=ctx.isoOligoCoverage(rows, F, "f_win");
  eq(c.nchg,0,"nchg"); eq(c.suggested,F,"unchanged"); eq(c.variants[0].pct,80,"pct"); });
T("coverage: fully conserved -> no variants, no change", ()=>{
  const c=ctx.isoOligoCoverage([mk(F),mk(F),mk(F)], F, "f_win");
  eq(c.variants.length,0,"variants"); eq(c.nchg,0,"nchg"); eq(c.mean,100,"mean"); });
T("coverage: <2 placed windows -> null", ()=>{ eq(ctx.isoOligoCoverage([mk(F)], F, "f_win"), null); });

// --- isoSummaryHtml: inclusivity / exclusivity breakdown ---
T("summary: counts + missed/cross-react lists, empty-safe", ()=>{
  const tg=[{role:"target",amplifies:true,probe_ident:100,probe_binds:true,acc:"A1"},
            {role:"target",amplifies:true,probe_ident:72,probe_binds:false,acc:"A2"},
            {role:"target",amplifies:false,acc:"A3"}];
  const ne=[{role:"neighbor",amplifies:false,acc:"N1"},
            {role:"neighbor",amplifies:true,probe_binds:true,acc:"N2"}];
  const h=ctx.isoSummaryHtml(tg,ne);
  has(h,"inclusivity"); has(h,"1/3","incl count"); has(h,"exclusivity"); has(h,"1/2","excl count");
  has(h,"A2","missed weak"); has(h,"A3","missed no-product"); has(h,"N2","cross-react");
  eq(ctx.isoSummaryHtml([],[]),"","empty-safe"); });

// --- isoLineageHtml: per-species detection table ---
T("lineage table: groups species, needs >=2 species", ()=>{
  const tg=[{role:"target",title:"Plasmodium relictum isolate 105 cytb",amplifies:true,probe_binds:true},
            {role:"target",title:"Plasmodium relictum isolate 88 cytb",amplifies:true,probe_binds:true},
            {role:"target",title:"Plasmodium vaughani isolate 427 cytb",amplifies:true,probe_binds:true}];
  const h=ctx.isoLineageHtml(tg); has(h,"relictum"); has(h,"vaughani"); has(h,"2/2","relictum detected");
  eq(ctx.isoLineageHtml([{role:"target",title:"Plasmodium relictum x"}]),"","single-species -> empty"); });

// --- isoCoverageHtml: reads iso_f/r/p, renders recommender ---
T("coverage html: forward recommender from iso_f", ()=>{
  el("iso_f").value="ACGTACGTAC"; el("iso_r").value=""; el("iso_p").value="";
  const rows=[mk(F),mk(F),mk(F),mk("ACGTATGTAC"),mk("ACGTATGTAC")];
  const h=ctx.isoCoverageHtml(rows); has(h,"coverage"); has(h,"ACGTAYGTAC","suggested forward");
  eq(ctx.isoCoverageHtml([mk(F)]),"","<2 -> empty"); });

// --- no-product reason categorizer + summary split by cause (v1.21.14 semantics) ---
T("reason categorizer + tag", ()=>{
  eq(ctx.isoReasonCat({reason:"forward region absent (best 60%)"}),"absent");
  eq(ctx.isoReasonCat({reason:"forward 3\u2032 mismatch 3 nt from end (best 95%)"}),"3prime");
  eq(ctx.isoReasonCat({reason:"primers bind but no convergent product in the size window"}),"size");
  eq(ctx.isoReasonTag({reason:"forward region absent"}),"region absent"); });
T("summary splits no-product by cause", ()=>{
  const tg=[{role:"target",amplifies:true,probe_ident:100,probe_binds:true,acc:"A1"},
            {role:"target",amplifies:false,acc:"A2",reason:"forward region absent (best 55%)"},
            {role:"target",amplifies:false,acc:"A3",reason:"reverse 3\u2032 mismatch 2 nt from end (best 95%)"}];
  const h=ctx.isoSummaryHtml(tg,[]);
  has(h,"region absent/partial"); has(h,"3\u2032 mismatch"); has(h,"A2"); has(h,"A3"); });

// --- results table filter + sort (v1.21.15) ---
T("isoSortFilter: filters", ()=>{
  const R=[{role:"target",acc:"T1",amplifies:true,probe_ident:100,probe_binds:true,product:150,amp_tm:89,f_ident:100,r_ident:100},
           {role:"target",acc:"T2",amplifies:true,probe_ident:60,probe_binds:false,product:150,amp_tm:90,f_ident:100,r_ident:100},
           {role:"target",acc:"T3",amplifies:false,f_ident:0,r_ident:100},
           {role:"neighbor",acc:"N1",amplifies:false,f_ident:0,r_ident:0},
           {role:"neighbor",acc:"N2",amplifies:true,probe_binds:true,product:150,amp_tm:88,f_ident:100,r_ident:100}];
  eq(ctx.isoSortFilter(R,{filter:"all",sort:"role",dir:1}).length,5);
  eq(ctx.isoSortFilter(R,{filter:"targets",sort:"role",dir:1}).length,3);
  eq(ctx.isoSortFilter(R,{filter:"neighbors",sort:"role",dir:1}).length,2);
  eq(ctx.isoSortFilter(R,{filter:"failures",sort:"role",dir:1}).map(x=>x.acc).join(","),"T2,T3"); // missed targets
  eq(ctx.isoSortFilter(R,{filter:"xreact",sort:"role",dir:1}).map(x=>x.acc).join(","),"N2"); });   // cross-reacts
T("isoSortFilter: sort by probe% asc/desc, nulls last", ()=>{
  const R=[{role:"target",acc:"A",probe_ident:90},{role:"target",acc:"B",probe_ident:60},{role:"target",acc:"C",probe_ident:null}];
  eq(ctx.isoSortFilter(R,{filter:"all",sort:"probe",dir:1}).map(x=>x.acc).join(","),"B,A,C");
  eq(ctx.isoSortFilter(R,{filter:"all",sort:"probe",dir:-1}).map(x=>x.acc).join(","),"A,B,C"); });

console.log("\n"+ok+" ok, "+bad.length+" failed"+(bad.length?": "+bad.join(", "):""));
process.exit(bad.length?1:0);
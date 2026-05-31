const fs=require("fs"),vm=require("vm");
const script=fs.readFileSync("static/index.html","utf8").match(/<script>([\s\S]*)<\/script>/)[1];
const VAL={std_in:"1e6, 20.1\n1e5, 23.4\n1e4, 26.7\n1e3, 30.0",email:"x@y.com"};
const store={}; const cache={};
function el(id){if(cache[id])return cache[id];cache[id]={_id:id,value:(VAL[id]!==undefined?VAL[id]:"ACGTACGTACGTACGT"),checked:false,innerHTML:"",textContent:"",style:{},classList:{add(){},remove(){},contains:()=>false},appendChild(){},addEventListener(){},click(){},querySelectorAll:()=>[],querySelector:()=>({value:"primer"})};return cache[id];}
const ctx={document:{getElementById:el,createElement:()=>el("c"+Math.random()),querySelectorAll:(s)=>(s==="nav a"||s==="section")?Array.from({length:16},(_,i)=>el("n"+i)):[]},window:{},console,Blob:function(){},URL:{createObjectURL:()=>"x"},fetch:async()=>({ok:true,status:200,json:async()=>({})}),localStorage:{getItem:k=>store[k]||null,setItem:(k,v)=>{store[k]=v;}},alert:()=>{},confirm:()=>true,setTimeout,JSON,Math,Object,Array,Date,parseInt,parseFloat,isFinite,Number,String,Boolean};
vm.createContext(ctx); vm.runInContext(script,ctx);
const RESP={
 "/api/pair":{f_tm:60,r_tm:60,pair_gap:0,fxr:-3,f_self:-4,r_self:-5,lint:[]},
 "/api/epcr":{forward_hits:1,reverse_hits:1,n_products:1,products:[{subject:"x",size:120,left:"F",right:"R",span:[10,130]}]},
 "/api/intron":{ok:true,verdict:"amplicon spans an exon-exon junction",info:"located on mRNA",amp_start:40,amp_end:176,amplicon:136,amp_located:true,junctions:[88,204],spanned:[88]},
 "/api/standard_curve":{efficiency_pct:99.2,slope:-3.33,r2:0.998,amp_factor:2,lod_practical:10,efficiency_ok:true,slope_ok:true,r2_ok:true,levels:[{quantity:"1e6",n:1,detected:1,mean_cq:20,sd_cq:null,detection_rate:100}],notes:"MIQE"},
 "/api/report":{html:"<html>r</html>",csv:"a,b\n1,2",n_assays:2},
 "/api/multiplex":{n_assays:2,n_oligos:4,threshold:-9,channel_conflicts:[{dye:"FAM",assays:["A","B"]}],cross_dimers:[{dg:-12,assay_a:"A",oligo_a:"F",assay_b:"B",oligo_b:"R"}],n_flagged:1},
};
ctx.api=async(path)=>RESP[path]||{};
(async()=>{let ok=0,bad=[];const T=async(label,fn)=>{try{await fn();ok++;console.log("  \u2713 "+label);}catch(e){bad.push(label);console.log("  \u2717 "+label+" -> "+e.message);}};
 await T("addAssay x2 + renderPanel",async()=>{ctx.addAssay({name:"FSJ IFNG",gene:"IFNG",organism:"Aphelocoma coerulescens",forward:"AGTCATTCTGATGTCGCTGATG",reverse:"ACCTGTCAGTGTTTTCAAGCA",probe:"TCATTTCTCTCTGTCCAGCCTGATAGCTTCTCT",amplicon:136,f_tm:60.1,r_tm:60.0,chemistry:"IDT PrimeTime (ZEN)"});ctx.addAssay({name:"Plas",gene:"Plasmodium cytb",forward:"TACCTGGACTWGTTTCATGG",reverse:"AAAGGATTTGTGCTACCTTG",probe:null,amplicon:157,chemistry:"low-Tm SYBR"});if(el("wk_out").innerHTML.indexOf("IFNG")<0)throw new Error("panel not rendered");});
 const html=el("wk_out").innerHTML; const id=(html.match(/setCurrent\('([^']+)'\)/)||[])[1];
 if(!id)throw new Error("no assay id parsed");
 await T("setCurrent + editAssay",async()=>{ctx.setCurrent(id);ctx.editAssay(id,"organism","Aphelocoma coerulescens");ctx.editAssay(id,"dye","FAM");});
 await T("toPair (prefill + run)",async()=>{ctx.toPair(id);if(el("p_f").value!=="AGTCATTCTGATGTCGCTGATG")throw new Error("p_f not prefilled");});
 await T("toEpcr (prefill)",async()=>{ctx.toEpcr(id);if(el("ep_f").value!=="AGTCATTCTGATGTCGCTGATG")throw new Error("ep_f not prefilled");});
 await T("toIntron (prefill)",async()=>{ctx.toIntron(id);if(el("i_gene").value!=="IFNG"||el("i_f").value!=="AGTCATTCTGATGTCGCTGATG")throw new Error("intron fields not prefilled");});
 await T("toConsAssay (prefill rows)",async()=>{ctx.toConsAssay(id);});
 await T("toOrderAssay (prefill rows)",async()=>{ctx.toOrderAssay(id);});
 await T("checkAssay (in-silico PCR)",async()=>{await ctx.checkAssay(id);});
 await T("doStd then attachStdToCurrent",async()=>{await ctx.doStd();ctx.attachStdToCurrent();if(el("wk_out").innerHTML.indexOf("val")<0)throw new Error("validation not attached/shown");});
 await T("genReport (download html+csv)",async()=>{await ctx.genReport();});
 await T("checkMultiplex (render)",async()=>{await ctx.checkMultiplex();if(el("wk_mx").innerHTML.indexOf("conflict")<0)throw new Error("multiplex not rendered");});
 await T("removeAssay",async()=>{ctx.removeAssay(id);});
 console.log("\n"+ok+" ok, "+bad.length+" failed"+(bad.length?": "+bad.join(", "):""));
})();

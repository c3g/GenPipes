var loadedobjects=""
var rootdomain="http://"+window.location.hostname

function ajaxpage(url, containerid){
	var page_request = false
	if (window.XMLHttpRequest) // if Mozilla, Safari etc
	page_request = new XMLHttpRequest()
	else if (window.ActiveXObject){ // if IE
	try {
		page_request = new ActiveXObject("Msxml2.XMLHTTP")
		} 
		catch (e){
			try{
				page_request = new ActiveXObject("Microsoft.XMLHTTP")
				}
				catch (e){}
				}
				}
				else
				return false
				page_request.onreadystatechange=function(){
					loadpage(page_request, containerid)
				}
					page_request.open('GET', url, true)
					page_request.send(null)
		}

function loadpage(page_request, containerid){
	if (page_request.readyState == 4 && (page_request.status==200 || window.location.href.indexOf("http")==-1))
	document.getElementById(containerid).innerHTML=page_request.responseText
	else if (page_request.readyState == 4) 
		document.getElementById(containerid).innerHTML="Loading...";
	else if (page_request.readyState == 3) 
		document.getElementById(containerid).innerHTML="In progress...";
		
	}
	

function toggleRowGroup(strVisibility, id){
	if(strVisibility == "visible"){
		
		if(document.all && !window.opera){ // && document.compatMode && document.compatMode == "CSS1Compat"){
			document.getElementById(id).style.display = "block";
		}else{
			document.getElementById(id).style.display = "table-row-group";
		}
	}else if(strVisibility == "collapse"){
		document.getElementById(id).style.display = "none";
	}
}
function toggleRowGroupVisibility(id){
	var list = document.getElementById("listFilter");
	//note first option is not considered (options[0]) since it's usually the default
	if(id == "all"){
		for(var i=1; i<list.options.length;i++){
			toggleRowGroup("visible", list.options[i].value);	
		}
	}else{
		//1. collapse all execept the requested one
		for(var i=1; i<list.options.length;i++){
			toggleRowGroup("collapse", list.options[i].value);	
		}
		//2. show requested row group
		toggleRowGroup("visible", id);	
	}
}

/********************************************
*	Toggle any element's display property 	*
*	to "none" (invisible)					*
*	@elNameRoot		is element id			*
********************************************/
function toggleVisibilityToNone(elNameRoot) {
	var el = document.getElementById(elNameRoot);
	var v = el.style.display = "none";
}

/********************************************
*	Toggle any element's display property 	*
*	to "block" (visible)					*
*	@elNameRoot		is element id			*
********************************************/
function toggleVisibilityToBlock(elNameRoot) {
	var el = document.getElementById(elNameRoot);
	var v = el.style.display = "block";
}
	
function toggleFAQVisibility(id){
	var list = document.getElementById("listFilter");
	//note first option is not considered (options[0]) since it's usually the default
	if(id == "all"){
		for(var i=1; i<list.options.length;i++){
			toggleVisibilityToBlock(list.options[i].value);	
		}
	}else{
		//1. collapse all execept the requested one
		for(var i=1; i<list.options.length;i++){
			toggleVisibilityToNone(list.options[i].value);	
		}
		//2. show requested row group
		toggleVisibilityToBlock(id);	
	}
}
function updateFilterList(id){
	var list = document.getElementById("listFilter");

	//note first option is not considered (options[0] = all or any) since it's usually the default
	for(var i=1; i<list.options.length;i++){
		if(list.options[i].value == id){
			list.options[i].selected = true;	
		}
	}
}

function toggleListBlockVisibility(id){
	var list = document.getElementsByTagName("li");

	if (list != null){
		for(var i=0; i<list.length;i++){
			if (list[i].getAttribute("name") == "mainListFilter"){
				toggleVisibilityToNone("div_" + list[i].id);	
			}
		}
	}
	toggleVisibilityToBlock("div_" + id);
}

function toggleItemCollectionVisibility(itemId){
	var el = document.getElementById(itemId);
	var list = document.getElementsByTagName("div");

	if (list != null){
		for(var i=0; i<list.length;i++){
			if(list[i].getAttribute("name") == el.getAttribute("name")){
				toggleVisibilityToNone(list[i].id);	
			}
		}
	}
	toggleVisibilityToBlock(itemId);	
}

function toggleItemCollectionToHidden(collectionName){
	var list = document.getElementsByTagName("div");
	if (list != null){
		for(var i=0; i<list.length;i++){
			if (list[i].getAttribute("name") == collectionName){
				toggleVisibilityToNone(list[i].id);	
			}
		}
	}
}

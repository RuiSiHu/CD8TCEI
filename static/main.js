String.prototype.trim=function()
{
	return this.replace(/(^\s*)(\s*$)/g, '');
}

String.prototype.ltrim=function()
{
	return this.replace(/(^\s*)/g,'');
}

String.prototype.rtrim=function()
{
	return this.replace(/(\s*$)/g,'');
}
function IgnoreSpaces(Str){
    	var ResultStr = "";
    	Temp=Str.split(" ");
    	for(i = 0; i < Temp.length; i++)
    		ResultStr +=Temp[i];
    	return ResultStr;
}

function clr(){
	form1.protein.value='';
	var obj=document.all.t2;
	obj.style.display="none";
	return false;
}

//清楚按钮方法
function clearVal(){
	//text="";
	//document.getElementById("textbyid").value=text;

	//$("#textbyid").val("");

	$("#textbyid").text("");
}

//清楚按钮方法二
function clearOne(){
	var buttonPlace=document.getElementsByClassName("enabled2");
	var textPlace=document.getElementsByClassName("textone");
	buttonPlace.onclick=function(){
		textPlace.value="";
	}

}

function checkSeqs()
{
	var seqstr=form1.protein.value;
 	seqstr = seqstr.replace(/\s/g,"");

	seqstr=IgnoreSpaces(seqstr);
 	//verify the protein of input is null
 	if(seqstr == null ||seqstr.length == 0){
 		alert('Please enter the protein sequence!');
 		document.form1.protein.focus()
 		return false;
 	}
	else{
		if(seqstr.length<50){
				alert("Your input sequence is less than 50aa and it probably is a fragment. Please input again!");
			form1.protein.focus();
			return(false);
		}
		else{
				str=seqstr.toLowerCase();
				var xnum=0;
				var amino="abcdefghiklmnpqrstvwy";
				for (var i=0; i<str.length;i++){
						var letter=str.charAt(i);
						if(letter=="x"){
								xnum=xnum+1;
						}
						else{
								if (amino.indexOf(letter) == -1){
										alert("Sorry,your sequence includes invalid character: " + letter +". Please see the example and input again!");
										form1.protein.focus();
										return(false);
								}
						}
				}
				if(xnum==str.length){
						alert("Sorry, you input a sequence with all X, we cannot go on working!");
						return(false);
				}
				if(xnum>4){
						question=confirm("Your input sequence include more than 4 X! Do you want to continue?");
						if(question==0){
								form1.protein.focus();
								return(false);
						}
				}
		}
	}
	form1.protein.value = seqstr;
	return(true);
}



function openwin(url)
{
    alert('hell0');
    alert(url);
    window.open(url, "newwindow", "height=500, width=600, top=0, left=0, toolbar=no, menubar=no, scrollbars=no, resizable=yes, location=no, status=no")

}



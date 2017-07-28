var nav_url = "https://raw.githubusercontent.com/sunericd/ISTools/master/ISWebsite/nav.json";

var tools_url = "https://raw.githubusercontent.com/sunericd/ISTools/master/ISWebsite/tools.json";

var dev_url="https://raw.githubusercontent.com/sunericd/ISTools/master/ISWebsite/developers.json";


$(document).ready(function(){

  /* Handlebars for navigation */
  $.getJSON(nav_url,
    function(data){
      var mysource = $('.nav').html();
      var mytemplate = Handlebars.compile(mysource);
      var myresult = mytemplate(data);
      $('.nav').html(myresult);

      /* Also distinguish current page */
      var mymainid = $('main').attr('id');
      var mynavid = '#navlink'+mymainid;
      $(mynavid).attr('id','iamhere');
    })

  /* Handlebars for tools on index.html*/
  $.getJSON(tools_url,
    function(data){
      var mysource = $('.tools_info').html();
      var mytemplate = Handlebars.compile(mysource);
      var myresult = mytemplate(data);
      $('.tools_info').html(myresult);
    })

    /* Handlebars for developers on About Us */
    $.getJSON(dev_url,
      function(data){
        var mysource = $('.person').html();
        var mytemplate = Handlebars.compile(mysource);
        var myresult = mytemplate(data);
        $('.person').html(myresult);
      })
    
})
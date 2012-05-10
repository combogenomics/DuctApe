kheader = '''
<head>
<style>
<!--
body {font-size:12px; font-family:verdana,arial,helvetica,sans-serif; background-color:#ffffff}
.text {font-size:12px; font-family:verdana,arial,helvetica,sans-serif}
h1 {font-size:24px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif; color:#005050}
h2 {font-size:18px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif; color:#005050}
h3 {font-size:16px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif}
h4 {font-size:14px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif}
h5 {font-size:12px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif; color:#444444}
b {font-size:12px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif}
th {font-size:12px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif}
td {font-size:12px; font-family:verdana,arial,helvetica,sans-serif}
input {font-size:12px; font-family:verdana,arial,helvetica,sans-serif}
select {font-size:12px; font-family:verdana,arial,helvetica,sans-serif}
pre {font-size:12px; font-family:courier,monospace}
a {text-decoration:none}
a:link {color:#003399}
a:visited {color:#003399}
a:hover {color:#33cc99}
font.title3 {color: #005050; font-size:14px; font-weight:bold; font-family:verdana,arial,helvetica,sans-serif}
hr.frame0 {border: 1px solid #f5e05f; color: #f5e05f}
div.poplay {
  position: absolute;
  padding: 2px;
  background-color: #ffff99;
  border-top: solid 1px #c0c0c0;
  border-left: solid 1px #c0c0c0;
  border-bottom: solid 1px #808080;
  border-right: solid 1px #808080;
  visibility: hidden;
}

span.popup
{
  font-weight: bold;
  color: #ffffff;
  white-space: nowrap;
}

form {
  margin: 0px;
}
-->
</style>
<script language="JavaScript">
var MSIE, Netscape, Opera, Safari, Firefox;

if(window.navigator.appName.indexOf("Internet Explorer") >= 0){
    MSIE = true;
}else if(window.navigator.appName == "Opera"){
    Opera = true;
}else if(window.navigator.userAgent.indexOf("Safari") >= 0){
    Safari = true;
}else if(window.navigator.userAgent.indexOf("Firefox") >= 0){
    Firefox = true;
    Netscape = true;
}else{
    Netscape = true;
}

function Component(id)
{
    this._component = document.getElementById(id);
    this._opacity_change_interval = 1;
    
    var opc = this._component.style.opacity;
    if(opc == "")
    {
        opc = 1;
    }
    this._opacity = opc * 100;
}


function _Component_ID()
{
    return this._component.id;
}
Component.prototype.id = _Component_ID;

function _Component_FontSize(size)
{
    if(_defined(size))
    {
        this._component.style.fontSize = size + "px";
    }
    else
    {
        return this._component.style.fontSize;
    }
}
Component.prototype.fontSize = _Component_FontSize;

function _Component_OpacityChangeInterval(interval)
{
    if(typeof(interval) == "undefined")
    {
        return this._opacity_change_interval;
    }
    else
    {
        this._opacity_change_interval = interval;
    }
}
Component.prototype.opacityChangeInterval = _Component_OpacityChangeInterval


function _Component_HTML(html)
{
    var component = this._component;
    
    if(typeof(html) == "undefined")
    {
        return component.innerHTML;
    }
    else
    {
        component.innerHTML = html;
    }
}

Component.prototype.HTML = _Component_HTML;

function _Component_BackgroundColor(color)
{
    this._component.style.backgroundColor = color;
}
Component.prototype.backgroundColor = _Component_BackgroundColor;

function _Component_BorderTop(border)
{
    if(_defined(border)){
        var comp = this._component;
        if(MSIE)
        {
            //comp.style.borderTop = border.color();
            //comp.style.border-top-style = border.style();
            //comp.style.border-top-width = border.width();
        }
        else
        {
            comp.style.borderTopColor = border.color();
            comp.style.borderTopStyle = border.style();
            comp.style.borderTopWidth = border.width() + "px";
        }
    }
}
Component.prototype.borderTop = _Component_BorderTop;

function _Component_BorderBottom(border)
{
    if(_defined(border)){
        var comp = this._component;
        if(MSIE)
        {
        }
        else
        {
            comp.style.borderBottomColor = border.color();
            comp.style.borderBottomStyle = border.style();
            comp.style.borderBottomWidth = border.width() + "px";
        }
    }
}
Component.prototype.borderBottom = _Component_BorderBottom;

function _Component_BorderLeft(border)
{
    if(_defined(border)){
        var comp = this._component;
        if(MSIE)
        {
        }
        else
        {
            comp.style.borderLeftColor = border.color();
            comp.style.borderLeftStyle = border.style();
            comp.style.borderLeftWidth = border.width() + "px";
        }
    }
}
Component.prototype.borderLeft = _Component_BorderLeft;

function _Component_BorderRight(border)
{
    if(_defined(border)){
        var comp = this._component;
        if(MSIE)
        {
        }
        else
        {
            comp.style.borderRightColor = border.color();
            comp.style.borderRightStyle = border.style();
            comp.style.borderRightWidth = border.width() + "px";
        }
    }
}
Component.prototype.borderRight = _Component_BorderRight;

function _Component_Border()
{
    var arg = _Component_Border.arguments;
    
    if(arg.length == 1)
    {
        this.borderTop(arg[0]);
        this.borderBottom(arg[0]);
        this.borderLeft(arg[0]);
        this.borderRight(arg[0]);
        
    }
    else if(arg.length == 2)
    {
        this.borderTop(arg[0]);
        this.borderBottom(arg[0]);
        this.borderLeft(arg[1]);
        this.borderRight(arg[1]);
        
    }else if(arg.length == 3)
    {
        this.borderTop(arg[0]);
        this.borderLeft(arg[1]);
        this.borderRight(arg[1]);
        this.borderBottom(arg[2]);
        
    }
    else if(arg.length == 4)
    {
        this.borderTop(arg[0]);
        this.borderRight(arg[1]);
        this.borderBottom(arg[2]);
        this.borderLeft(arg[3]);
    }
}
Component.prototype.border = _Component_Border;

function _Component_X(x)
{
    var component = this._component;
    
    if(typeof(x) == "undefined")
    {
        var ret = (MSIE) ? component.style.pixelLeft : parseInt(component.style.left);
        return ret;
    }
    else
    {
        if(MSIE)
        {
            component.style.pixelLeft = x;
        }
        else if(Opera)
        {
            component.style.left = x;
        }
        else
        {
            component.style.left = x + "px";
        }
    }
}
Component.prototype.x = _Component_X;

function _Component_Y(y)
{
    var component = this._component;
    
    if(typeof(y) == "undefined")
    {
        var ret = (MSIE) ? component.style.pixelTop : parseInt(component.style.top);
        return ret;
    }else
    {
        if(MSIE)
        {
            component.style.pixelTop = y;
        }
        else if(Opera)
        {
            component.style.top = y;
        }
        else
        {
            component.style.top = y + "px";
        }
    }
}
Component.prototype.y = _Component_Y;

function _Component_Move(x, y)
{
    this.x(x);
    this.y(y);
}
Component.prototype.move = _Component_Move;

function _Component_Width(width)
{
    var component = this._component;
    
    if(typeof(width) == "undefined")
    {
        var ret = (MSIE) ? component.style.pixelWidth : parseInt(component.style.width);
        return ret;
    }
    else
    {
        if(MSIE)
        {
            component.style.pixelWidth = width;
        }
        else if(Opera)
        {
            component.style.width = width;
        }
        else
        {
            component.style.width = width + "px";
        }
    }
}
Component.prototype.width = _Component_Width;

function _Component_Height(height)
{
    var component = this._component;
    
    if(typeof(height) == "undefined")
    {
        var ret = (MSIE) ? component.style.pixelWidth : parseInt(component.style.width);
        return ret;
    }
    else
    {
        if(MSIE)
        {
            component.style.pixelHeight = height;
        }
        else if(Opera)
        {
            component.style.height = height;
        }
        else
        {
            component.style.height = height + "px";
        }
    }
}
Component.prototype.height = _Component_Height;

function _Component_Size(width, height)
{
    this.width(width);
    this.height(height);
}
Component.prototype.size = _Component_Size;

function _Component_Visible(visible)
{
    var component = this._component;
    
    if(typeof(visible) == "undefined")
    {
        return (component.style.visibility == "visible") ? true : false;
    }
    else
    {
        if(MSIE || Safari || Firefox || Opera)
        {
            if(visible)
            {
                component.style.visibility = "visible";
                this._opacityStep = 10;
                this._opacity = 0;
            }
            else
            {
                this._opacityStep = -10;
                this._opacity = this.opacity();
            }
            
            _addComponent(this);
            
            this.changeOpacity();
        }
        else
        {
            component.style.visibility = (visible) ? "visible" : "hidden";
        }
    }
}
Component.prototype.visible = _Component_Visible;

function _Component_ChangeOpacity()
{
    var opacity = this._opacity + this._opacityStep;
    
    this.opacity(opacity);
    
    if(opacity >= 100)
    {
        return;
    }
    else if(opacity <= 0)
    {
        this._component.style.visibility = "hidden";
        return
    }
    else
    {
        var interval = this._opacity_change_interval;
        setTimeout("_triggerChangeOpacity('" + this.id() + "')", interval);
    }
}
Component.prototype.changeOpacity = _Component_ChangeOpacity;

function _Component_Opacity(opacity)
{
    if(typeof(opacity) == "undefined")
    {
        return this._opacity;
    }
    else
    {
        this._opacity = opacity;
        
        var component = this._component;
        component.style.opacity = opacity / 100;
        component.style.mozOpacity = opacity / 100;
        component.style.filter = "alpha(opacity=" + opacity + ")";
    }
}
Component.prototype.opacity = _Component_Opacity;

var _component_list = new Array();

function _addComponent(component)
{
    var id = component.id();
    _component_list[id] = component;
}

function _triggerChangeOpacity(id)
{
    var component = _component_list[id];
    component.changeOpacity();
}

function _defined(val)
{
    return (typeof(val) != "undefined") ? true : false;
}

function Border()
{
    this._width = 1;
    this._style = "solid";
    this._color = "#000000";
}

function _Border_Color(color)
{
    if(!_defined(color)){
        return this._color;
    }else{
        this._color = color;
    }
}
Border.prototype.color = _Border_Color;

function _Border_Style(style)
{
    if(!_defined(style)){
        return this._style;
    }else{
        this._style = style;
    }
}
Border.prototype.style = _Border_Style;

function _Border_Width(width)
{
    if(!_defined(width)){
        return this._width;
    }else{
        this._width = width;
    }
}
Border.prototype.width = _Border_Width;

document.onmousemove = _documentMouseMove;

var _mousePosX = 0;
var _mousePosY = 0;

function _documentMouseMove(evt)
{
    _mousePosX = _getEventX(evt);
    _mousePosY = _getEventY(evt);
}

function _getEventX(evt)
{
    var ret;
    if(Netscape){
        ret = evt.pageX;
    }else if(MSIE){
        ret = event.x + getPageXOffset();
    }else if(Safari){
        ret = event.x + getPageXOffset();
    }else{
        ret = evt.x;
    }

    return ret;
}

function _getEventY(evt)
{
    var ret;

    if(Netscape){
        ret = evt.pageY;
    }else if(MSIE){
        ret = event.y + getPageYOffset();
    }else if(Safari){
        ret = event.y + getPageYOffset();
    }else{
        ret = event.y;
    }

    return ret;
}

function getCurrentMouseX()
{
    return _mousePosX;
}

function getCurrentMouseY()
{
    return _mousePosY
}

function getPageXOffset()
{
    var ret;
    if(Safari || Opera){
        ret = document.body.scrollLeft;
    }else{
        if(document.body.scrollLeft > 0){
            ret = document.body.scrollLeft;
        }else{
            ret = document.documentElement.scrollLeft;
        }
    }

    return ret;
}

function getPageYOffset()
{
    var ret;
    if(Safari || Opera){
        ret = document.body.scrollTop;
    }else{
        if(document.body.scrollTop > 0){
            ret = document.body.scrollTop;
        }else{
            ret = document.documentElement.scrollTop;
        }
    }

    return ret;
}            
            
            
var timer = 0;
var p_entry, p_title, p_bgcolor;
function popupTimer(entry, title, bgcolor)
{
  p_entry = entry;
  p_title = title;
  p_bgcolor = bgcolor;

  if(timer == 0){
    var func = "showThumbnail()";
    timer = setTimeout(func, 1200);
  }
}


function showThumbnail()
{

  var url = "";
  if(p_entry.match(/^[A-Z]+\d+$/))
  {
    url = "http://www.genome.jp/kegg/misc/thumbnail/" + p_entry + ".gif";
  }
  else if(p_entry.match(/(\d+)$/))
  {
    url = "http://www.genome.jp/kegg/misc/thumbnail/map" + RegExp.$1 + ".gif";
  }

  var html = "";

  html += '<img src="' + url + '" alt="Loading...">';

  var x = getCurrentMouseX();
  var y = getCurrentMouseY();

  var layer = new Component("poplay");
  layer.backgroundColor(p_bgcolor);
  layer.HTML(html);
  layer.move(x, y+40);
  layer.visible(true);

  timer = 0;
}


function hideMapTn(){
  var layer = new Component("poplay");
  layer.visible(false);

  if(timer != 0){
    clearTimeout(timer);
    timer = 0;
  }
}
</script>
</head>
'''


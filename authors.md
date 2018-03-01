---
layout: page
title: About the Authors
permalink: /authors/
---

<div id="authors">
{% for author in site.data.authors %}
<h5 id="{{ username }}">{{ author[1].name }}</h5>
<p align="center">
<img src='"{{ username }}">{{ author[1].assets}}'
/>
</p>
<p id="{{ username }}">{{ author[1].location }}</p>
<p id="{{ username }}">{{ author[1].bio }}</p>
<p More info about this author at: href="{{ username }}">{{ author[1].url }}</p>
<ul class="posts">
{% for post in site.posts %}
{% if author[1].username == post.author %}
{% if post.title != null %}
<li itemscope><span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}" itemprop="datePublished">{{ post.date | date: "%B %d, %Y" }}</time></span> &raquo; <a href="{{ site.baseurl }}{{ post.url | remove: '/'}}">{{ post.title }}</a></li>
{% endif %}
{% endif %}
{% endfor %}
</ul>
<hr>
{% endfor %}
</div>

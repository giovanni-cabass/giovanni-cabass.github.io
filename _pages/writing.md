---
layout: page
title: Notes
permalink: /writing/
---

{% if site.posts.size > 0 %}
<ul class="post-list">
  {% for post in site.posts %}
  <li class="post-entry">
    <div class="post-title"><a href="{{ post.url | relative_url }}">{{ post.title }}</a></div>
    <div class="post-date">{{ post.date | date: "%B %-d, %Y" }}</div>
    {% if post.excerpt %}<div class="post-excerpt">{{ post.excerpt }}</div>{% endif %}
  </li>
  {% endfor %}
</ul>
{% else %}
<p class="writing-empty">No posts yet — check back soon.</p>
{% endif %}

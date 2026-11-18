# Publish the project website with GitHub Pages

The website is static HTML and CSS in `docs/index.html` and `docs/site.css`.
The simulation video is reused from `docs/media/`. No website dependencies,
build tools, or server-side application are required.

## Preview locally

From the repository root:

```sh
python3 -m http.server 8000 --directory docs
```

Open `http://localhost:8000/`. Stop the preview with Ctrl+C.

## Enable GitHub Pages

After the website files are committed and pushed:

1. Open the repository's **Settings → Pages** on GitHub.
2. Under **Build and deployment**, choose **Deploy from a branch**.
3. Select **master**, the **/docs** folder, and **Save**.
4. Wait for GitHub's Pages deployment to complete. The settings page will show
   the published URL.

The expected project URL is:

`https://saptarshi7889.github.io/Conservative-Finite-Difference-Method-for-GPE-of-1d-Bose-Polaron/`

This URL becomes live only after Pages is enabled and deployment succeeds.
Subsequent pushes to the publishing source update the site automatically.
`docs/.nojekyll` keeps the HTML/CSS site independent of Jekyll processing.

Keep the version and release-download links in `docs/index.html` in sync when
publishing a new solver release. Hosting the website does not run numerical
simulations in the visitor's browser; visitors download and build the solver.

[GitHub's publishing-source documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).
